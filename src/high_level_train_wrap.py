import optuna
from catboost import CatBoostRegressor
from src.utils import eval_metrics

from sklearn.decomposition import PCA

import numpy as np

from src.multitask_nn import RegressionModel
from torch.utils.data import TensorDataset

from argparse import Namespace

from lightning.pytorch.callbacks import Callback, ModelCheckpoint
from lightning.pytorch.loggers import TensorBoardLogger
from lightning import seed_everything, Trainer
from torch.utils.data import DataLoader
import torch

from datetime import datetime

def objective(trial, X_train, X_test, y_train, y_test):
    params = {
        "iterations": trial.suggest_int("iterations", 300, 1000),
        "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3),
        "depth": trial.suggest_int("depth", 4, 10),
        "l2_leaf_reg": trial.suggest_float("l2_leaf_reg", 1e-2, 10.0),
        "bagging_temperature": trial.suggest_float("bagging_temperature", 0.0, 1.0),
        "random_strength": trial.suggest_float("random_strength", 1e-3, 5.0),
        "border_count": trial.suggest_int("border_count", 32, 255),
        "verbose": 0,
        "loss_function": "RMSE",
        "task_type": "CPU"
    }

    model = CatBoostRegressor(**params)
    model.fit(X_train, y_train, eval_set=(X_test, y_test), early_stopping_rounds=30)

    y_pred = model.predict(X_test)
    metrics = eval_metrics(y_test, y_pred, "regression")
    trial.set_user_attr("RMSE", metrics["RMSE"])  # Store RMSE inside the trial
    return metrics["R2"]

def best_trial_callback(study, trial):
    if study.best_trial.number == trial.number:
        rmse = trial.user_attrs.get("RMSE")
        print(f"🏆 New Best Trial {trial.number}: R2={trial.value:.5f}, RMSE={rmse:.5f}")

def train_optuna_catboost(X_train, y_train, X_val, y_val, n_trials: int = 50, n_jobs: int = 8) -> dict:
    study = optuna.create_study(direction="maximize")
    study.optimize(
        lambda trial: objective(trial, X_train, X_val, y_train, y_val),
        n_trials=n_trials,
        n_jobs=n_jobs,
        show_progress_bar=True,
        callbacks=[best_trial_callback]
    )

    # Train final model with best params
    best_params = study.best_params
    best_params["loss_function"] = "RMSE"
    best_params["verbose"] = 0
    
    return best_params

def train_regression_nn(X_train: np.ndarray, y_train: np.ndarray, X_val: np.ndarray, y_val: np.ndarray, cfg: Namespace):
    X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train, dtype=torch.float32)

    X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
    y_val_tensor = torch.tensor(y_val, dtype=torch.float32)
    
    train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
    val_dataset = TensorDataset(X_val_tensor, y_val_tensor)

    # persistent_workers=True reduces overhead of creating workers
    train_loader = DataLoader(train_dataset, batch_size=cfg.batch_size, shuffle=True, num_workers=cfg.num_workers, persistent_workers=cfg.persistent_workers)
    val_loader = DataLoader(val_dataset, batch_size=cfg.batch_size, shuffle=False, num_workers=cfg.num_workers, persistent_workers=cfg.persistent_workers)
    
    input_dim = X_train_tensor.shape[1]
    model = RegressionModel(
        input_dim=input_dim,

        hidden_dim=cfg.hid_dim,
        depth=cfg.depth,
        use_residual=cfg.use_residual,
        drop=cfg.drop,

        lr=cfg.lr,
        scheduler_gamma=cfg.scheduler_gamma,
    )

    current_date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    logger = TensorBoardLogger(
        save_dir="tb_logs/", name=cfg.project_name
    )

    checkpoint_callback = ModelCheckpoint(
        dirpath=f"checkpoints/{cfg.project_name}-{current_date}",
        filename="{epoch:02d}-{val_loss:.4f}",
        save_top_k=1,
        monitor="V_tot",
        mode="min",
        save_last=True,
    )

    trainer = Trainer(
        logger=logger,
        callbacks=[checkpoint_callback],
        max_epochs=cfg.max_epochs,
        gradient_clip_val=1.0,
    )
    trainer.fit(model, train_loader, val_loader)
    
    return model