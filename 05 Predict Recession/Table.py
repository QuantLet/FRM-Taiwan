import streamlit as st
import pandas as pd
from itertools import product
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.stats import rankdata

st.set_page_config(layout="wide")

df = pd.read_csv('Performance.csv')

config_order = ['Sector', 'Window', 'Selection', 'Regression', 'Scale', 'Response', 'Variable', 'Target', 'Goal', 'Method', 'Feature']


st.title('Table Configurations')

col1, col2 = st.columns([0.25, 0.75], vertical_alignment="bottom")
with col1:
    st.subheader('Columns:')
with col2:
    col_vars = st.multiselect(
        '↓  Select variables for columns.',
        options=config_order,
        default=['Sector', 'Window']
    )
    col_vars = sorted(col_vars, key=lambda x: config_order.index(x))

if col_vars:
    max_columns = 4
    num_vars = len(col_vars)
    for i in range(0, num_vars, max_columns):
        cols = st.columns(min(max_columns, max(num_vars, max_columns)))
        for j, var in enumerate(col_vars[i:i + max_columns]):
            values = [value[0].upper() + value[1:] if isinstance(value, str) else value for value in sorted(df[var].unique())]
            values = [str(value) for value in values]
            with cols[j]:
                st.markdown(f"**{var} ({len(values)}):** {' / '.join(values)}")

change_combinations = list(product(*(df[var].unique() for var in col_vars)))

col1, col2 = st.columns([0.25, 0.75], vertical_alignment="bottom")
with col1:
    st.subheader('Rows:')
with col2:
    row_vars = st.multiselect(
        '↓  Select variables for rows.',
        options=config_order,
        default=['Variable']
    )
    row_vars = sorted(row_vars, key=lambda x: config_order.index(x))

if row_vars:
    max_columns = 4
    num_vars = len(row_vars)
    for i in range(0, num_vars, max_columns):
        cols = st.columns(min(max_columns, max(num_vars, max_columns)))
        for j, var in enumerate(row_vars[i:i + max_columns]):
            values = [value[0].upper() + value[1:] if isinstance(value, str) else value for value in sorted(df[var].unique())]
            values = [str(value) for value in values]
            with cols[j]:
                st.markdown(f"**{var} ({len(values)}):** {' / '.join(values)}")

change_combinations = list(product(*(df[var].unique() for var in row_vars)))


col1, col2 = st.columns([0.25, 0.75], vertical_alignment="bottom")
with col1:
    st.subheader('Scores:')
with col2:
    score_order = ['In-AUC', 'In-Brier', 'In-LL', 'Out-AUC', 'Out-Brier', 'Out-LL']
    score_vars = st.multiselect(
        '↓  Select multiple performance metrics.',
        options=score_order,
        default=['In-AUC', 'Out-AUC']
    )
    score_vars = sorted(score_vars, key=lambda x: score_order.index(x))

st.markdown('---')

if row_vars and col_vars:
    row_combinations = list(product(*(df[var].unique() for var in row_vars)))
    col_combinations = list(product(*(df[var].unique() for var in col_vars)))

    for layer in range(len(col_vars)):
        header_cols = st.columns(len(col_combinations) + len(row_vars))
        if layer == 0:
            for r in range(len(row_vars)):
                with header_cols[r]:
                    st.markdown(f"**{row_vars[r]}**")
        for i, col_combo in enumerate(col_combinations):
            formatted_value = col_combo[layer][0].upper() + col_combo[layer][1:] if isinstance(col_combo[layer], str) else col_combo[layer]
            with header_cols[i + len(row_vars)]:
                st.markdown(f"**{formatted_value}**")

    for row_combo in row_combinations:
        row_cols = st.columns(len(col_combinations) + len(row_vars))
        for layer in range(len(row_vars)):
            with row_cols[layer]:
                formatted_value = row_combo[layer][0].upper() + row_combo[layer][1:] if isinstance(row_combo[layer], str) else row_combo[layer]
                st.markdown(f"**{formatted_value}**")
        for j, col_combo in enumerate(col_combinations):
            filtered_df = df.copy()
            for var, value in zip(row_vars + col_vars, row_combo + col_combo):
                filtered_df = filtered_df[filtered_df[var] == value]

            avg_values = filtered_df[score_vars].mean()
            with row_cols[j + len(row_vars)]:
                st.markdown(f"{avg_values.mean():.3f}")


