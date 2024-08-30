import streamlit as st
import pandas as pd
from itertools import product
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.stats import rankdata

st.set_page_config(layout="wide")

df = pd.read_csv('Performance.csv')

st.title('Performance')

col1, col2 = st.columns([0.25, 0.75], vertical_alignment="bottom")
with col1:
    st.subheader('Configurations')
with col2:
    config_order = ['Sector', 'Window', 'Selection', 'Regression', 'Scale', 'Response', 'Variable', 'Target', 'Goal', 'Method', 'Feature']
    change_vars = st.multiselect(
        '↓  Choose multiple variables for comparison.',
        options=config_order,
        default=['Sector', 'Window']
    )
    change_vars = sorted(change_vars, key=lambda x: config_order.index(x))

if change_vars:
    max_columns = 4
    num_vars = len(change_vars)
    for i in range(0, num_vars, max_columns):
        cols = st.columns(min(max_columns, max(num_vars, max_columns)))
        for j, var in enumerate(change_vars[i:i + max_columns]):
            values = [value[0].upper() + value[1:] if isinstance(value, str) else value for value in sorted(df[var].unique())]
            values = [str(value) for value in values]
            with cols[j]:
                st.markdown(f"**{var} ({len(values)}):** {' / '.join(values)}")

change_combinations = list(product(*(df[var].unique() for var in change_vars)))

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

col1, col2 = st.columns([0.25, 0.75], vertical_alignment="bottom")
with col1:
    st.subheader('Statistics:')
with col2:
    stat_order = ['Box Plot', 'Avg', 'Med', 'Std', 'Max', 'Min']
    stat_vars = st.multiselect(
        '↓  Choose statistical metrics to display.',
        options=stat_order,
        default=['Box Plot', 'Avg']
    )
stat_vars = sorted(stat_vars, key=lambda x: stat_order.index(x))

st.markdown('---')

results = {}
for combo in change_combinations:
    filtered_df = df.copy()
    for var, value in zip(change_vars, combo):
        filtered_df = filtered_df[filtered_df[var] == value]
    
    result = {}
    for score_var in score_vars:
        result[score_var] = {}
        if 'Avg' in stat_vars:
            result[score_var]['Avg'] = filtered_df[score_var].mean()
        if 'Med' in stat_vars:
            result[score_var]['Med'] = filtered_df[score_var].median()
        if 'Std' in stat_vars:
            result[score_var]['Std'] = filtered_df[score_var].std()
        if 'Max' in stat_vars:
            result[score_var]['Max'] = filtered_df[score_var].max()
        if 'Min' in stat_vars:
            result[score_var]['Min'] = filtered_df[score_var].min()
    
    results[combo] = result

ranks = {}
for score_var in score_vars:
    for stat_var in stat_vars:
        if stat_var == 'Box Plot':
            continue
        if score_var in ['In-AUC', 'Out-AUC']:
            ascending = False
        elif score_var in ['In-Brier', 'In-LL', 'Out-Brier', 'Out-LL']:
            ascending = True
        else:
            ascending = False

        if stat_var in ['Std', 'Min']:
            ascending = True

        stat_values = [results[combo][score_var][stat_var] for combo in change_combinations]
        ranked_values = rankdata(stat_values, method='min' if ascending else 'dense')

        if not ascending:
            ranked_values = len(ranked_values) - ranked_values + 1

        for combo, rank in zip(change_combinations, ranked_values):
            if combo not in ranks:
                ranks[combo] = {}
            if score_var not in ranks[combo]:
                ranks[combo][score_var] = {}
            ranks[combo][score_var][stat_var] = rank

def md_head(score_var):
    return f'''
    <div style="display: inline-block; text-align: left;">
        <div style="font-size: 20px;margin-bottom: 4px;">
            <span style="font-weight: 800;">{score_var}</span>
        </div>
    </div><br>
    '''

def md_body(stat_var, score_var, stat_value, rank, total):
    midpoint = total / 2
    color = "red" if rank <= midpoint else "green"
    
    return f'''
    <div style="display: inline-block; text-align: left;">
        <div style="font-size: 16px; font-weight: 600;">
            <span style="font-weight: 800;">{stat_var}.</span>
        </div>
        <div style="font-size: 32px; font-weight: 200; color: dimgray;margin-top: -12px;">
            {stat_value:.3f}
        </div>
        <div style="font-size: 14px; font-weight: 800; color: {color}; margin-top: -12px;margin-bottom: 10px;">
            # {rank}
        </div>
    </div><br>
    '''


if len(change_vars) > 0 and len(score_vars) > 0:
    total_combinations = len(change_combinations)
    for i in range(0, total_combinations, 6):
        columns = st.columns(7)
        with columns[0]:
            for var in change_vars:
                st.markdown(f"##### {var}")
        
        for j, combo in enumerate(change_combinations[i:i + 7]):
            if j < len(columns) - 1:
                with columns[j + 1]:
                    for value in combo:
                        if isinstance(value, str):
                            value = value[0].upper() + value[1:]
                        st.markdown(f"##### {value}")
                    for score_var in score_vars:
                        sub_cols = st.columns(2)
                        with sub_cols[0]:
                            metric_html = ""
                            metric_html += md_head(score_var)
                            for stat_var in stat_vars:
                                if stat_var == 'Box Plot':
                                    continue
                                stat_value = results[combo][score_var][stat_var]
                                rank = ranks[combo][score_var][stat_var]
                                metric_html += md_body(stat_var, score_var, stat_value, rank, total_combinations)
                            st.markdown(metric_html, unsafe_allow_html=True)
                        
                        if 'Box Plot' in stat_vars:
                            with sub_cols[1]:
                                fig, ax = plt.subplots(figsize=(4,20)) 
                                sns.boxplot(data=df[df[change_vars[0]] == combo[0]][score_var], ax=ax,
                                            color='gray', width=0.4,
                                            linewidth=8, linecolor='black',
                                            flierprops={'marker':'o', 'markersize':10},
                                            boxprops={'linewidth':4})
                                if score_var in ['In-AUC', 'Out-AUC']:
                                    ax.set_ylim(-0.05, 1.05)
                                elif score_var in ['In-Brier', 'Out-Brier']:
                                    ax.set_ylim(0, df[score_var].max())
                                elif score_var in ['In-LL', 'Out-LL']:
                                    ax.set_ylim(-0.5,df[score_var].max())
                                ax.set_ylabel('')
                                ax.spines['top'].set_linewidth(5)
                                ax.spines['right'].set_linewidth(5)
                                ax.spines['left'].set_linewidth(5)
                                ax.spines['bottom'].set_linewidth(5)
                                ax.grid(True, axis='y', linestyle='--', linewidth=6, color='gray', alpha=0.6)
                                font_prop = FontProperties(weight='bold', size=36)
                                for label in ax.get_yticklabels():
                                    label.set_fontproperties(font_prop)
                                plt.tight_layout()
                                st.pyplot(fig)

                        st.write('')
                            
        st.markdown('---')
else:
    st.warning("Please select at least one configuration/score.")