import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import io
import plotly.express as px

@st.cache_data(show_spinner=False)
def prepare_h5ad(file_bytes):
    adata = sc.read_h5ad(io.BytesIO(file_bytes))
    processed_adata = adata.copy()
    
    sc.pp.normalize_total(
        processed_adata,
        target_sum=10_000,
    )
    sc.pp.log1p(processed_adata)
    
    sc.pp.highly_variable_genes(
        processed_adata,
        n_top_genes=min(2000, processed_adata.n_vars),
    )
    
    number_of_hvgs = processed_adata.var["highly_variable"].sum()
    
    number_of_pcs = min(
        50,
        processed_adata.n_obs - 1,
        int(number_of_hvgs) - 1,
    )
    
    sc.tl.pca(
        processed_adata,
        n_comps=number_of_pcs,
        use_highly_variable=True,
    )
    
    sc.pp.neighbors(
        processed_adata,
        n_neighbors=15,
        n_pcs=number_of_pcs,
    )

    try:
        sc.tl.leiden(processed_adata, key_added="clusters")
    except ImportError:
        processed_adata.obs["clusters"] = "Clustering dependency unavailable"
    
    sc.tl.umap(processed_adata)
    
    return adata, processed_adata

st.title("Cellalyze: Single-Cell Gene Expression Explorer")

st.write(
    "This app will help explore single-cell gene-expression data, "
    "including clusters, marker genes, and UMAP visualizations."
)

st.sidebar.header("Controls")

uploaded_file = st.sidebar.file_uploader(
    "Upload a dataset",
    type=["h5ad", "csv", "gz"],
    accept_multiple_files=False,
)

if uploaded_file is None:
    st.warning("Please upload a dataset to get started.")
else:
    st.success(f"Uploaded file: {uploaded_file.name}")

    preview_rows = st.sidebar.slider(
        "Number of preview rows",
        min_value=5,
        max_value=20,
        value=5,
        step=1,
    )

    show_column_names = st.sidebar.checkbox("Show column names", value=True)

    filename = uploaded_file.name.lower()

    if filename.endswith(".h5ad"):
        st.write("This is an AnnData file.")
        
        with st.spinner("Preparing dataset..."):
            adata, processed_adata = prepare_h5ad(
                uploaded_file.getvalue()
            )
            
        st.success("Loaded successfully")
        sample = adata.X[
            : min(500, adata.n_obs),
            : min(500, adata.n_vars),
        ]

        if hasattr(sample, "toarray"):
            sample = sample.toarray()

        st.write("Expression matrix type:", type(adata.X).__name__)
        st.write("Raw matrix available:", adata.raw is not None)
        st.write("Sample minimum:", float(np.min(sample)))
        st.write("Sample maximum:", float(np.max(sample)))
        st.write(
            "Sample contains only whole numbers:",
            bool(np.allclose(sample, np.round(sample))),
        )
        
        processed_sample = processed_adata.X[
            : min(500, processed_adata.n_obs),
            : min(500, processed_adata.n_vars),
        ]
        
        if hasattr(processed_sample, "toarray"):
            processed_sample = processed_sample.toarray()
            
        st.write(
            "Processed sample maximum:",
            float(np.max(processed_sample)),
        )
        st.write(
            "Processed sample contains only whole numbers:",
            bool(
                np.allclose(
                    processed_sample,
                    np.round(processed_sample),
                )
            ),
        )
        
        st.success("Normalization and log transformation completed.")
        
        number_of_hvgs = processed_adata.var["highly_variable"].sum()
        
        st.write(
            "Number of highly variable genes:",
            int(number_of_hvgs),
        )
        
        st.write(
            "PCA coordinates shape:",
            processed_adata.obsm["X_pca"].shape,
        )
        
        st.write(
            "Neighbor graph created:",
            "neighbors" in processed_adata.uns,
        )
        
        st.write(
            "UMAP coordinates shape:",
            processed_adata.obsm["X_umap"].shape,
        )

        color_mode = st.sidebar.radio(
            "Color UMAP by",
            ["Gene expression", "Cell metadata"],
        )

        if color_mode == "Gene expression":
            plot_color = st.sidebar.selectbox(
                "Choose a gene",
                processed_adata.var_names.tolist(),
            )
        else:
            metadata_columns = processed_adata.obs.columns.tolist()
            if not metadata_columns:
                st.info("This dataset has no cell metadata columns.")
                plot_color = None
            else:
                plot_color = st.sidebar.selectbox(
                    "Choose a metadata column",
                    metadata_columns,
                    index=(
                        metadata_columns.index("clusters")
                        if "clusters" in metadata_columns
                        else 0
                    ),
                )

        umap_coordinates = processed_adata.obsm["X_umap"]
        
        umap_df = pd.DataFrame(
            {
                "UMAP 1": umap_coordinates[:, 0],
                "UMAP 2": umap_coordinates[:, 1],
            },
            index=processed_adata.obs_names,
        )
        
        if color_mode == "Gene expression":
            gene_expression = processed_adata[
                :, plot_color
            ].X
            
            if hasattr(gene_expression, "toarray"):
                gene_expression = gene_expression.toarray().ravel()

            umap_df[plot_color] = gene_expression
        elif plot_color is not None:
            umap_df[plot_color] = processed_adata.obs[plot_color].to_numpy()
                
        st.subheader("UMAP Data Preview")
        st.dataframe(umap_df.head())
        
        interactive_fig = px.scatter(
            umap_df,
            x="UMAP 1",
            y="UMAP 2",
            color=plot_color,
            hover_name=umap_df.index,
            title=f"Interactive UMAP by {plot_color or 'cell'}",
        )
        
        interactive_fig.update_traces(
            marker={"size": 5, "opacity": 0.7}
        )
        
        st.plotly_chart(
            interactive_fig,
            use_container_width=True,
        )
        
        st.write("Example gene names:", adata.var_names[:10].tolist())
        st.write("Available embeddings:", list(processed_adata.obsm.keys()))
        st.write("Cell metadata columns:", list(processed_adata.obs.columns))
        st.write("Available layers:", list(adata.layers.keys()))
        st.write("Analysis keys:", list(processed_adata.uns.keys()))
        st.write(f"Number of cells: {adata.n_obs}")
        st.write(f"Number of genes: {adata.n_vars}")

        if plot_color is not None:
            fig = sc.pl.umap(
                processed_adata,
                color=plot_color,
                show=False,
                return_fig=True,
            )

            st.pyplot(fig)
            
    elif filename.endswith(".csv.gz"):
        st.write("This is a gzipped CSV file.")
        data = pd.read_csv(uploaded_file, compression="gzip", comment="#")

        st.subheader("Data Preview")
        st.dataframe(data.head(preview_rows))

        st.write(
            f"The dataset contains {data.shape[0]} rows "
            f"and {data.shape[1]} columns."
        )

        if show_column_names:
            st.subheader("Column Names")
            st.write(data.columns.tolist())

        selected_column = st.sidebar.selectbox("Choose a column", data.columns)

        st.write(f"You selected: {selected_column}")

        st.write("Values in the selected column:")
        st.write(data[selected_column].head())

        st.subheader("Simple Plot")
        if pd.api.types.is_numeric_dtype(data[selected_column]):
            fig, ax = plt.subplots()
            ax.hist(
                data[selected_column].dropna(),
                bins=20,
                color="steelblue",
                edgecolor="black",
            )
            ax.set_title(f"Distribution of {selected_column}")
            ax.set_xlabel(selected_column)
            ax.set_ylabel("Count")
            st.pyplot(fig)
        else:
            counts = data[selected_column].value_counts().head(20)
            fig, ax = plt.subplots()
            counts.plot(kind="bar", ax=ax, color="skyblue")
            ax.set_title(f"Frequency of {selected_column}")
            ax.set_xlabel(selected_column)
            ax.set_ylabel("Count")
            st.pyplot(fig)

    elif filename.endswith(".csv"):
        st.write("This is a CSV file.")
        data = pd.read_csv(uploaded_file, comment="#")
        st.subheader("Data Preview")
        st.dataframe(data.head(preview_rows))

        st.write(
            f"The dataset contains {data.shape[0]} rows "
            f"and {data.shape[1]} columns."
        )

        if show_column_names:
            st.subheader("Column Names")
            st.write(data.columns.tolist())

        selected_column = st.sidebar.selectbox("Choose a column", data.columns)

        st.write(f"You selected: {selected_column}")

        st.write("Values in the selected column:")
        st.write(data[selected_column].head())

        st.subheader("Simple Plot")
        if pd.api.types.is_numeric_dtype(data[selected_column]):
            fig, ax = plt.subplots()
            ax.hist(
                data[selected_column].dropna(),
                bins=20,
                color="steelblue",
                edgecolor="black",
            )
            ax.set_title(f"Distribution of {selected_column}")
            ax.set_xlabel(selected_column)
            ax.set_ylabel("Count")
            st.pyplot(fig)
        else:
            counts = data[selected_column].value_counts().head(20)
            fig, ax = plt.subplots()
            counts.plot(kind="bar", ax=ax, color="skyblue")
            ax.set_title(f"Frequency of {selected_column}")
            ax.set_xlabel(selected_column)
            ax.set_ylabel("Count")
            st.pyplot(fig)

    else:
        st.warning(
            "Unsupported file type. Please upload a .h5ad, .csv, or .csv.gz file."
        )


