import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def plot_T_from_distance_for_nn(model, val_loader, threshold=0.5, figsize=(8, 5)):
    result_logits, result_reg = [], []

    for batch in tqdm(val_loader, desc="Computing T(distance) plot"):
        class_logits, reg_output = model(batch[0])
        result_logits.append(class_logits.detach().cpu().numpy())
        result_reg.append(reg_output.detach().cpu().numpy())

    result_logits = np.concatenate(result_logits).flatten()
    result_reg = np.concatenate(result_reg).flatten()

    # Apply sigmoid to get probabilities
    probs = 1 / (1 + np.exp(-result_logits))
    # Get predicted classes (1 if prob > 0.5, else 0)
    pred_classes = (probs > threshold).astype(int)

    # Plot with different colors for each class
    plt.figure(figsize=figsize)
    
    # Colors and linear regression for each class
    for name, cls, color in [('Liquid', 0, 'blue'), ('Solid', 1, 'red')]:
        mask = pred_classes == cls
        
        # Plot scatter points
        plt.scatter(
            result_logits[mask], result_reg[mask], 
            c=color, 
            label=name,
            alpha=0.7
        )
        
        # Perform linear regression for this class
        if np.sum(mask) > 1:  # Need at least 2 points for regression
            x = result_logits[mask]
            y = result_reg[mask]
            slope, intercept = np.polyfit(x, y, 1)
            
            # Plot regression line
            x_range = np.array([min(x), max(x)])
            plt.plot(x_range, slope * x_range + intercept, 
                     color=color, linestyle='--', 
                     label=f'{name} fit: {slope:.2f}x + {intercept:.2f}')
            
    plt.ylabel('Temperature, C')
    plt.xlabel('Distance to separating plane')
    plt.grid(alpha=0.3)
    
    plt.legend()
    plt.show()
