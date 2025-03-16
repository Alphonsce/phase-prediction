from typing import Tuple, List

import torch
import torch.nn as nn
import torch.nn.functional as F

from sklearn.model_selection import train_test_split

from lightning import LightningModule, Trainer

from torchmetrics.classification import Accuracy, F1Score, AUROC
from torchmetrics import R2Score, MeanSquaredError, MeanAbsoluteError

from torch.utils.data import TensorDataset


class MultiTaskModel(LightningModule):
    def __init__(
        self,
        input_dim,
        hidden_dim=256,
        drop=0.2,
        lr=1e-3,

        cl_loss_coef=1,
        reg_loss_coef=1,
    ):
        super().__init__()

        self.lr = lr
        self.cl_loss_coef = cl_loss_coef
        self.reg_loss_coef = reg_loss_coef

        self.save_hyperparameters()

        self.shared = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(drop),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(drop),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU()
        )


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

    def training_step(self, batch, batch_idx):
        x, labels, temp = batch

        logits, predictions = self(x)

        classification_loss = F.binary_cross_entropy_with_logits(logits, labels.unsqueeze(1).float())
        regression_loss = F.mse_loss(predictions, temp.unsqueeze(1).float())

        loss = self.cl_loss_coef * classification_loss + self.reg_loss_coef * regression_loss

        self.log("T_tot", loss, prog_bar=True)
        self.log("T_cl", classification_loss, prog_bar=True)
        self.log("T_reg", regression_loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, labels, temp = batch

        logits, predictions = self(x)

        classification_loss = F.binary_cross_entropy_with_logits(logits, labels.unsqueeze(1).float())
        regression_loss = F.mse_loss(predictions, temp.unsqueeze(1).float())

        loss = self.cl_loss_coef * classification_loss + self.reg_loss_coef * regression_loss

        self.log("V_tot", loss, prog_bar=True)

        probs = torch.sigmoid(logits).squeeze()
        preds = (probs > 0.5).long()
        labels_long = labels.long()

        self.val_accuracy(preds, labels_long)
        self.val_f1(preds, labels_long)
        self.val_roc_auc(probs, labels_long)
        self.val_r2_class(probs, labels.float())

        self.val_mse(predictions.squeeze(), temp)
        self.val_mae(predictions.squeeze(), temp)
        self.val_r2_reg(predictions.squeeze(), temp)

        return loss

    def on_validation_epoch_end(self):
        # Classification:
        self.log("metrics/classification/acc", self.val_accuracy.compute(), prog_bar=True)
        self.log("metrics/classification/f1_macro", self.val_f1.compute())
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

def create_datasets(X, labels, temp, train_size=0.8) -> Tuple[TensorDataset]:
    dataset_components = train_test_split(X, labels, temp, train_size=train_size)

    # np.ndarray -> torch.tensor:
    X_train, X_val, y_train, y_val, temp_train, temp_val = list(map(lambda x: torch.tensor(x, dtype=torch.float32), dataset_components))
    
    train_data = TensorDataset(X_train, y_train, temp_train)
    val_data = TensorDataset(X_val, y_val, temp_val)

    return train_data, val_data