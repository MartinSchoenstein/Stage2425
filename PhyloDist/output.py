def to_csv(distance_matrix, p_value=None, name="matrix.csv", pvalue_name="p_value.csv"):
    
    distance_matrix.to_csv(name, index=True)
    if p_value is not None:
        p_value.to_csv(pvalue_name, index=True)
        return distance_matrix, p_value
    else:
        return distance_matrix
