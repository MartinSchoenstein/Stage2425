import pandas as pd

def to_csv(distance_matrix, index, columns):
    
    distance_matrix = pd.DataFrame(distance_matrix, index=index, columns=columns)
    distance_matrix.to_csv(f"{distance_matrix.name}.csv", index=True)
    
