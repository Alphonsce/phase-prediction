from typing import Tuple, List

import torch
import torch.nn as nn
import torch.nn.functional as F

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from lightning import LightningModule, Trainer

from torchmetrics.classification import Accuracy, F1Score, AUROC
from torchmetrics import R2Score, MeanSquaredError, MeanAbsoluteError

from torch.utils.data import TensorDataset

class MLPBlock(nn.Module):
    def __init__(self, input_dim, output_dim, drop=0.2):
        super().__init__()
        self.block = nn.ModuleList([
            nn.Linear(input_dim, output_dim),
            nn.BatchNorm1d(output_dim),
            nn.ReLU(),
            nn.Dropout(drop)
        ])
        
    def forward(self, x):
        for layer in self.block:
            x = layer(x)
        return x

class ResidualBlock(nn.Module):
    def __init__(self, block: nn.ModuleList):
        super().__init__()
        self.block = block
        
    def forward(self, x):
        return x + self.forward_block(x)
        
    def forward_block(self, x):
        for layer in self.block:
            x = layer(x)
        return x

class MultiTaskModel(LightningModule):
    def __init__(
        self,
        input_dim,
        hidden_dim=256,
        depth=3,
        use_residual=False,
        drop=0.2,
        lr=1e-3,
        cl_loss_coef=1,
        reg_loss_coef=1,
    ):
        super().__init__()

        self.lr = lr
        self.cl_loss_coef = cl_loss_coef
        self.reg_loss_coef = reg_loss_coef
        self.use_residual = use_residual

        self.save_hyperparameters()

        # Build shared network with configurable depth
        shared_layers = []
        shared_layers.append(MLPBlock(input_dim, hidden_dim, drop))
        
        for _ in range(depth - 1):
            if use_residual:
                block = nn.ModuleList([
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.BatchNorm1d(hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(drop)
                ])
                shared_layers.append(ResidualBlock(block))
            else:
                shared_layers.append(MLPBlock(hidden_dim, hidden_dim, drop))
            
        self.shared = nn.Sequential(*shared_layers)

        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

        self.val_accuracy = Accuracy(task='binary')
        self.val_f1 = F1Score(task='binary')
        self.val_roc_auc = AUROC(task='binary')
        self.val_r2_class = R2Score()

        self.val_mse = MeanSquaredError()
        self.val_mae = MeanAbsoluteError()
        self.val_r2_reg = R2Score()

    def forward(self, x):
        shared_out = self.shared(x)
        class_logits = self.classifier(shared_out)
        reg_output = self.regressor(shared_out)

        return class_logits, reg_output
    
    def compute_masked_regression_loss(self, temp_pred, temp_true):
        '''
        Computes the regression loss only for samples with not NaN Temperatures.
        '''
        mask = ~torch.isnan(temp_true)
        if mask.sum() > 0:
            valid_pred_temp = temp_pred[mask.unsqueeze(1)].squeeze()
            valid_temp = temp_true[mask]
            return F.mse_loss(valid_pred_temp, valid_temp)
        return torch.tensor(0.0, device=self.device)

    def training_step(self, batch, batch_idx):
        x, labels, temp_true = batch

        logits, temp_pred = self(x)

        classification_loss = F.binary_cross_entropy_with_logits(logits, labels.unsqueeze(1).float())
        regression_loss = self.compute_masked_regression_loss(temp_pred, temp_true)

        loss = self.cl_loss_coef * classification_loss + self.reg_loss_coef * regression_loss

        self.log("T_tot", loss, prog_bar=True)
        self.log("T_cl", classification_loss, prog_bar=True)
        self.log("T_reg", regression_loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, labels, temp_true = batch

        logits, temp_pred = self(x)

        classification_loss = F.binary_cross_entropy_with_logits(logits, labels.unsqueeze(1).float())
        regression_loss = self.compute_masked_regression_loss(temp_pred, temp_true)

        loss = self.cl_loss_coef * classification_loss + self.reg_loss_coef * regression_loss

        self.log("V_tot", loss, prog_bar=True)

        probs = torch.sigmoid(logits).squeeze()
        preds = (probs > 0.5).long()
        labels_long = labels.long()

        self.val_accuracy(preds, labels_long)
        self.val_f1(preds, labels_long)
        self.val_roc_auc(probs, labels_long)
        self.val_r2_class(probs, labels.float())

        mask = ~torch.isnan(temp_true)
        if mask.sum() > 0:
            self.val_mse(temp_pred.squeeze()[mask], temp_true[mask])
            self.val_mae(temp_pred.squeeze()[mask], temp_true[mask])
            self.val_r2_reg(temp_pred.squeeze()[mask], temp_true[mask])
        else:
            self.log("metrics/regression/mse", float('nan'))
            self.log("metrics/regression/mae", float('nan'))
            self.log("metrics/regression/r2_reg", float('nan'))

        return loss

    def on_validation_epoch_end(self):
        # Classification:
        self.log("metrics/classification/acc", self.val_accuracy.compute())
        self.log("metrics/classification/f1", self.val_f1.compute(), prog_bar=True)
        self.log("metrics/classification/roc_auc", self.val_roc_auc.compute())
        self.log("metrics/classification/r2_class", self.val_r2_class.compute())
        # Regression:
        self.log("metrics/regression/mse", self.val_mse.compute())
        self.log("metrics/regression/mae", self.val_mae.compute(), prog_bar=True)
        self.log("metrics/regression/r2_reg", self.val_r2_reg.compute())

        self.val_accuracy.reset()
        self.val_f1.reset()
        self.val_roc_auc.reset()
        self.val_r2_class.reset()
        
        self.val_mse.reset()
        self.val_mae.reset()
        self.val_r2_reg.reset()

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        return optimizer

def create_datasets(X, labels, temp, train_size=0.8, use_norm=False) -> Tuple[TensorDataset]:
    X_train, X_val, y_train, y_val, temp_train, temp_val = train_test_split(X, labels, temp, train_size=train_size)
    
    if use_norm:
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_val = scaler.transform(X_val)
    
    X_train, X_val, y_train, y_val, temp_train, temp_val = list(map(lambda x: torch.tensor(x, dtype=torch.float32), [X_train, X_val, y_train, y_val, temp_train, temp_val]))
    
    train_data = TensorDataset(X_train, y_train, temp_train)
    val_data = TensorDataset(X_val, y_val, temp_val)

    return train_data, val_data