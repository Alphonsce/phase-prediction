from typing import Tuple, List

import torch
import torch.nn as nn
import torch.nn.functional as F

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from lightning import LightningModule, Trainer

from torchmetrics.classification import Accuracy, F1Score, AUROC
from torchmetrics import R2Score, MeanSquaredError, MeanAbsoluteError
from math import log

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
        scheduler_gamma=1,
        cl_loss_coef=1,
        pos_label_weight=100,
        reg_loss_coef=1,
        use_learnable_loss_weights=False,  # New argument for learnable weights
    ):
        super().__init__()

        # Store hyperparameters
        self.lr = lr
        self.scheduler_gamma = scheduler_gamma
        self.cl_loss_coef = cl_loss_coef
        self.reg_loss_coef = reg_loss_coef
        self.use_residual = use_residual
        self.pos_label_weight = pos_label_weight
        self.use_learnable_loss_weights = use_learnable_loss_weights

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

        # Classification head
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

        # Regression head
        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

        # Define learnable parameters for loss weighting if enabled
        if self.use_learnable_loss_weights:
            self.log_var_cl = nn.Parameter(torch.tensor(log(cl_loss_coef)))  # Log variance for classification
            self.log_var_reg = nn.Parameter(torch.tensor(log(reg_loss_coef)))  # Log variance for regression

        # Metrics
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
        Computes the regression loss only for samples with non-NaN temperatures.
        '''
        mask = ~torch.isnan(temp_true)
        if mask.sum() > 0:
            valid_pred_temp = temp_pred[mask.unsqueeze(1)].squeeze()
            valid_temp = temp_true[mask]
            return F.mse_loss(valid_pred_temp, valid_temp.squeeze())
        return torch.tensor(0.0, device=self.device)

    def training_step(self, batch, batch_idx):
        x, labels, temp_true = batch
        logits, temp_pred = self(x)

        # Compute individual losses
        classification_loss = F.binary_cross_entropy_with_logits(
            logits, labels.unsqueeze(1).float(), pos_weight=torch.tensor([self.pos_label_weight])
        )
        regression_loss = self.compute_masked_regression_loss(temp_pred, temp_true)

        # Compute total loss based on weighting method
        if self.use_learnable_loss_weights:
            total_loss = (
                (torch.exp(-self.log_var_cl) * classification_loss + self.log_var_cl) +
                (torch.exp(-self.log_var_reg) * regression_loss + self.log_var_reg)
            )
            self.log("T_weight_cl", torch.exp(-self.log_var_cl), prog_bar=False)
            self.log("T_weight_reg", torch.exp(-self.log_var_reg), prog_bar=False)
        else:
            total_loss = self.cl_loss_coef * classification_loss + self.reg_loss_coef * regression_loss

        # Logging
        self.log("T_tot", total_loss, prog_bar=True)
        self.log("T_cl", classification_loss, prog_bar=True)
        self.log("T_reg", regression_loss, prog_bar=True)
        return total_loss

    def validation_step(self, batch, batch_idx):
        x, labels, temp_true = batch
        logits, temp_pred = self(x)

        # Compute individual losses
        classification_loss = F.binary_cross_entropy_with_logits(
            logits, labels.unsqueeze(1).float(), pos_weight=torch.tensor([self.pos_label_weight])
        )
        regression_loss = self.compute_masked_regression_loss(temp_pred, temp_true)

        # Compute total loss based on weighting method
        if self.use_learnable_loss_weights:
            weight_cl = torch.exp(-self.log_var_cl)
            weight_reg = torch.exp(-self.log_var_reg)
            total_loss = weight_cl * classification_loss + weight_reg * regression_loss
            self.log("V_weight_cl", weight_cl, prog_bar=False)
            self.log("V_weight_reg", weight_reg, prog_bar=False)
        else:
            total_loss = self.cl_loss_coef * classification_loss + self.reg_loss_coef * regression_loss

        # Logging
        self.log("V_tot", total_loss, prog_bar=True)
        self.log("V_cl", classification_loss, prog_bar=False)
        self.log("V_reg", regression_loss, prog_bar=False)

        # Classification metrics
        probs = torch.sigmoid(logits).squeeze()
        preds = (probs > 0.5).long()
        labels_long = labels.long()
        self.val_accuracy(preds, labels_long)
        self.val_f1(preds, labels_long)
        self.val_roc_auc(probs, labels_long)
        self.val_r2_class(probs, labels.float())

        # Regression metrics
        mask = ~torch.isnan(temp_true)
        if mask.sum() > 0:
            self.val_mse(temp_pred.squeeze()[mask], temp_true[mask])
            self.val_mae(temp_pred.squeeze()[mask], temp_true[mask])
            self.val_r2_reg(temp_pred.squeeze()[mask], temp_true[mask])
        else:
            self.log("metrics/regression/mse", float('nan'))
            self.log("metrics/regression/mae", float('nan'))
            self.log("metrics/regression/r2_reg", float('nan'))

        return total_loss

    def on_validation_epoch_end(self):
        # Classification metrics
        self.log("metrics/classification/acc", self.val_accuracy.compute())
        self.log("metrics/classification/f1", self.val_f1.compute(), prog_bar=True)
        self.log("metrics/classification/roc_auc", self.val_roc_auc.compute())
        self.log("metrics/classification/r2_class", self.val_r2_class.compute())
        # Regression metrics
        mse = self.val_mse.compute()
        self.log("metrics/regression/mse", mse)
        self.log("metrics/regression/mae", self.val_mae.compute())
        self.log("metrics/regression/rmse", torch.sqrt(mse), prog_bar=True)
        self.log("metrics/regression/r2_reg", self.val_r2_reg.compute(), prog_bar=True)

        # Reset metrics
        self.val_accuracy.reset()
        self.val_f1.reset()
        self.val_roc_auc.reset()
        self.val_r2_class.reset()
        self.val_mse.reset()
        self.val_mae.reset()
        self.val_r2_reg.reset()

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=self.scheduler_gamma)
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": scheduler,
                "interval": "epoch",
                "frequency": 1
            }
        }

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