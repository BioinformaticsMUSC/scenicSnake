#!/usr/bin/env python3
"""
Generate comprehensive HTML report for SCENIC analysis
"""

import re
import pandas as pd
import numpy as np
import json
import os
from datetime import datetime
import plotly.graph_objects as go
import plotly
import scanpy as sc

def main():
    # Get parameters from snakemake
    auc_file = snakemake.input.auc_matrix
    regulons_file = snakemake.input.regulons
    rss_file = snakemake.input.rss_scores
    plot_files = snakemake.input.plots
    output_file = snakemake.output[0]
    metadata_file = snakemake.input.metadata

    cell_type_column = snakemake.params.cell_type_column

    print("Generating SCENIC analysis report...")
    
    # Load data
    auc_df = pd.read_csv(auc_file, index_col=0)
    
    regulons = get_regs(regulon_csv_file=regulons_file)
    
    rss_df = pd.read_csv(rss_file, index_col=0)

    adata = sc.read_h5ad(metadata_file)
    # split_column = snakemake.config['loom_preparation']['split_condition']
    # adata = adata[adata.obs[split_column] == snakemake.params.split_value].copy()
    
    # Ensure cell order matches
    common_cells = list(set(auc_df.index) & set(adata.obs_names))
    auc_df = auc_df.loc[common_cells, :]
    adata = adata[common_cells]
    auc_df_ct = pd.concat([auc_df, adata.obs[cell_type_column]], axis=1)
    
    # Generate HTML report
    html_content = generate_html_report(auc_df, auc_df_ct, regulons, rss_df, plot_files, cell_type_column)

    # Save report
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"Report saved to {output_file}")

def get_regs(regulon_csv_file):
    """convert regulon csv file to regulon dict/json format"""
    df = pd.read_csv(regulon_csv_file,
    index_col=0, 
	skiprows=2, 
	names=['MotifID', 'AUC', 'NES','MotifSimilarityQvalue','OrthologousIdentity','Annotation', 'Context', 'TargetGenes', 'RankAtMax'], 
	header=0)
    reg_json_data = {
        f"{row.name}(+)": [gene for gene, _ in parse_TargetGenes(row['TargetGenes'])]
        for _, row in df.iterrows()
    }
    return reg_json_data

def parse_TargetGenes(x):
	reg_patt = r"\('(\w*)',\s(?:np.float64\()?(\d*.\d*)(?:\))?\)"
	genes = re.findall(reg_patt, x)
	return genes

def generate_html_report(auc_df, auc_df_ct, regulons, rss_df, plot_files, cell_type_column):
    """Generate HTML report content"""
    
    # Calculate summary statistics
    n_regulons = len(regulons)
    n_cells = auc_df.shape[0]
    mean_targets_per_regulon = np.mean([len(genes) for genes in regulons.values()])
    
    # Find top regulons by activity and specificity
    mean_activity = auc_df.mean(axis=0).sort_values(ascending=False)
    top_active_regulons = mean_activity.head(10)
    
    max_rss = rss_df.max(axis=0).sort_values(ascending=False)
    top_specific_regulons = max_rss.head(10)
    
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SCENIC Analysis Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
        }}
        
        .section {{
            background: white;
            padding: 25px;
            margin-bottom: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        
        .stat-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            border-left: 4px solid #667eea;
        }}
        
        .stat-number {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        
        .stat-label {{
            color: #666;
            font-size: 0.9em;
            margin-top: 5px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        
        th {{
            background-color: #f8f9fa;
            font-weight: 600;
        }}
        
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        
        .regulon-list {{
            columns: 2;
            column-gap: 30px;
        }}
        
        .regulon-item {{
            break-inside: avoid;
            margin-bottom: 10px;
            padding: 10px;
            background: #f8f9fa;
            border-radius: 5px;
        }}
        
        .timestamp {{
            color: #666;
            font-size: 0.9em;
            text-align: center;
            margin-top: 30px;
        }}
        
        h1, h2, h3 {{
            color: #333;
        }}
        
        h2 {{
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        
        .alert {{
            padding: 15px;
            margin-bottom: 20px;
            border: 1px solid transparent;
            border-radius: 4px;
        }}
        
        .alert-info {{
            color: #31708f;
            background-color: #d9edf7;
            border-color: #bce8f1;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 SCENIC Analysis Report</h1>
        <p>Single-Cell Regulatory Network Inference and Clustering</p>
    </div>
    
    <div class="section">
        <h2>📊 Analysis Summary</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-number">{n_regulons}</div>
                <div class="stat-label">Regulons Discovered</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{n_cells:,}</div>
                <div class="stat-label">Cells Analyzed</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{mean_targets_per_regulon:.1f}</div>
                <div class="stat-label">Mean Targets per Regulon</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{len(rss_df.index)}</div>
                <div class="stat-label">Cell Types</div>
            </div>
        </div>
        
        <div class="alert alert-info">
            <strong>Analysis Complete!</strong> SCENIC identified {n_regulons} gene regulatory networks 
            across {n_cells:,} cells. The regulons contain an average of {mean_targets_per_regulon:.1f} 
            target genes each.
        </div>
    </div>
    
    <div class="section">
        <h2>🔥 Top Active Regulons</h2>
        <p>Regulons with highest average activity across all cells:</p>
        <table>
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Regulon</th>
                    <th>Mean AUC Score</th>
                    <th>Target Genes</th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Add top active regulons
    for i, (regulon, score) in enumerate(top_active_regulons.items(), 1):
        target_genes = ", ".join(regulons.get(regulon, {}))
        html += f"""
                <tr>
                    <td>{i}</td>
                    <td><strong>{regulon}</strong></td>
                    <td>{score:.3f}</td>
                    <td>{target_genes}</td>
                </tr>
"""
    
    html += f"""
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>🎯 Most Cell-Type Specific Regulons</h2>
        <p>Regulons with highest specificity to particular cell types:</p>
        <table>
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Regulon</th>
                    <th>Max RSS Score</th>
                    <th>Most Specific Cell Type</th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Add top specific regulons
    for i, (regulon, score) in enumerate(top_specific_regulons.items(), 1):
        if regulon in rss_df.T.index:
            most_specific_cell_type = rss_df.T.loc[regulon].idxmax()
        else:
            most_specific_cell_type = 'N/A'
        
        html += f"""
                <tr>
                    <td>{i}</td>
                    <td><strong>{regulon}</strong></td>
                    <td>{score:.2f}</td>
                    <td>{most_specific_cell_type}</td>
                </tr>
"""
    html += f"""
            </tbody>
        </table>
    </div>
"""
    html += f"""
    <div class="section">
        <h2>🎯 Top Regulons Specific to each Cell type</h2>
"""
    for celltype in rss_df.index:
        top10_regs = rss_df.loc[celltype].nlargest(10)
        html += f"""
        <h4>{celltype}:</h4>
        <table>
            <thead>
                <tr>
                <th>Rank</th>
                <th>Regulon</th>
                <th>RSS Score</th>
                </tr>
            </thead>
            <tbody>
"""
        for i, (regulon, score) in enumerate(top10_regs.items(), 1):
            html += f"""
                <tr>
                    <td>{i}</td>
                    <td><strong>{regulon}</strong></td>
                    <td>{score:.2f}</td>
                </tr>
            """
        html += """
                </tbody>
                </table>
"""
    html += f"""
    </div>
    
    <div class="section">
        <h2>📈 Visualization Gallery</h2>
        <p>Interactive plots showing regulon activity patterns and cell type specificity.</p>
"""
    
    # Add plot sections
    plot_titles = {
        'regulon_heatmap.pdf': 'Regulon Activity Heatmap',
        'umap_regulon_activity.pdf': 'UMAP with Regulon Activity',
        'rss_plot.pdf': 'Regulon Specificity Scores'
    }
    ##create plotly heatmap
    reg_dict = {}
    for ct in rss_df.index:
        reg_dict[ct] = list(rss_df.T[ct].sort_values(ascending=False)[:5].index)
    auc_means = auc_df_ct.groupby(cell_type_column).mean().T
    top_regs_list = []
    for c in reversed(auc_means.columns):
        top_cc = reg_dict[c]
        top_cc_sorted = auc_means.loc[top_cc,c].sort_values(ascending=True).index
        for r in top_cc_sorted:
            if r not in top_regs_list:
                top_regs_list.append(r)
    hm_data = auc_means.loc[top_regs_list,:]
    fig = go.Figure(data=go.Heatmap(
        z=hm_data.values,
        x=hm_data.columns,
        y=hm_data.index,
        colorscale='Reds',
        colorbar=dict(title='AUC Score')
    ))
    fig.update_layout(title='Regulon Activity Heatmap', xaxis_title='Cell Types', yaxis_title='Regulons')
    
    #plotly.offline.plot(fig, filename='results/plots/regulon_activity_heatmap.html')
    html += plotly.io.to_html(fig, full_html=False, include_plotlyjs='cdn')
    for plot_file in plot_files:
        plot_name = os.path.basename(plot_file)
        plot_title = plot_titles.get(plot_name, plot_name)
        
        html += f"""
        <div class="plot-container">
            <h3>{plot_title}</h3>
            <p><em>Plot saved as: {plot_name}</em></p>
        </div>
"""
    
    html += f"""
    </div>
    
    <div class="section">
        <h2>📋 All Discovered Regulons</h2>
        <p>Complete list of {n_regulons} regulons identified in this analysis:</p>
        <div class="regulon-list">
"""
    
    # Add all regulons
    for regulon_name, regulon_data in sorted(regulons.items()):
        n_targets = len(regulon_data)
        html += f"""
            <div class="regulon-item">
                <strong>{regulon_name}</strong><br>
                <small>{n_targets} target genes: {", ".join(regulon_data) if n_targets > 0 else "None"}</small>
            </div>
"""
    
    html += f"""
        </div>
    </div>
    
    <div class="section">
        <h2>📁 Output Files</h2>
        <p>The following files were generated during this analysis:</p>
        <ul>
            <li><strong>AUC Matrix:</strong> <code>results/scenic/auc_matrix.csv</code></li>
            <li><strong>Binary Activity:</strong> <code>results/scenic/binary_regulon_activity.csv</code></li>
            <li><strong>Regulons:</strong> <code>results/scenic/regulons.json</code></li>
            <li><strong>RSS Scores:</strong> <code>results/scenic/rss_scores.csv</code></li>
            <li><strong>Visualization Plots:</strong> <code>results/plots/</code></li>
        </ul>
    </div>
    
    <div class="section">
        <h2>🔬 Methods</h2>
        <p><strong>SCENIC</strong> (Single-Cell rEgulatory Network Inference and Clustering) is a computational method 
        to infer gene regulatory networks from single-cell RNA-seq data and to identify cell states.</p>
        
        <h3>Analysis Steps:</h3>
        <ol>
            <li><strong>Co-expression Network:</strong> GENIE3 algorithm identifies co-expressed gene modules</li>
            <li><strong>Regulon Discovery:</strong> cisTarget performs motif enrichment to identify transcription factor regulons</li>
            <li><strong>Activity Scoring:</strong> AUCell calculates regulon activity scores for each cell</li>
            <li><strong>Cell State Identification:</strong> Regulon activity patterns define cellular states</li>
        </ol>
        
        <h3>Citations:</h3>
        <ul>
            <li>Aibar et al. (2017). SCENIC: single-cell regulatory network inference and clustering. <em>Nature Methods</em></li>
            <li>Van de Sande et al. (2020). A scalable SCENIC workflow for single-cell gene regulatory network analysis. <em>Nature Protocols</em></li>
        </ul>
    </div>
    
    <div class="timestamp">
        <p>Report generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>SCENIC Snakemake Workflow v1.0</p>
    </div>
</body>
</html>
"""
    
    return html

if __name__ == "__main__":
    main()
