import pandas as pd
import numpy as np
import plotly.graph_objects as go

def load_data(filename):
    df = pd.read_csv(filename, delim_whitespace=True, names=['x', 'y', 'z'])
    return df


def function_surface(x, y):
    return np.sin(np.pi*x) * np.sin(np.pi*y)
    

def plot_3d(data):

    x_range = np.linspace(data['x'].min(), data['x'].max(), 50)
    y_range = np.linspace(data['y'].min(), data['y'].max(), 50)
    X, Y = np.meshgrid(x_range, y_range)
    Z = function_surface(X, Y)
    
    scatter = go.Scatter3d(
        x=data['x'], y=data['y'], z=data['z'],
        mode='markers',
        marker=dict(size=5, color=data['z'], colorscale='Viridis', opacity=1)
    )
    

    surface = go.Surface(
        x=X, y=Y, z=Z,
        colorscale=[[0, 'rgb(255, 204, 255)'], [1, 'rgb(255, 204, 255)']],
        opacity=0.5,
        showscale=False
    )
    
    fig = go.Figure(data=[scatter, surface])
    
    fig.update_layout(
        title='Интерактивный 3D-график',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        )
    )
    
    fig.show()


if __name__ == "__main__":
    filename = '/Users/mihail/Desktop/somecode9.1/somecode9.1/solution.txt'  # Замените на имя вашего файла
    data = load_data(filename)
    plot_3d(data)
