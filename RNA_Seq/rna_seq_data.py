import warnings
from types import NoneType
import os
from os import PathLike
from pathlib import Path
import hashlib
import json
from datetime import datetime
from importlib.metadata import version as pkg_version

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import seaborn as sns  # noqa: F401  (registers the "rocket_r" colormap)
from adjustText import adjust_text

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

class RNASeq_Data():
    def __init__(self, filepath: PathLike, sep: str, 
                 index_col_name: str, cache_dir: str | None,
                 sample_number_cat_name: NoneType | str, sample_number_cats: NoneType | dict
        ):
        # Data sources
        # If filepath is false, it's not loaded yet.
        self.data_sources = {filepath: False}
        self.index_col_name = index_col_name
        self.df = self._load_csv(filepath, sep, sample_number_cat_name, sample_number_cats)

        # Helper variables for setting up and maintaining comparisons
        self.currently_comparing = False
        self.current_comparison = None
        self.comparison_table = None
        self.dds = None
        self.ds = None
        self.comparison_p = None
        self.comparison_results = None
        self.comparison_from_cache = False
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None

    def _load_csv(self, filepath, sep, sample_number_cat_name: NoneType | str, sample_number_cats: NoneType | dict, drop_objects: bool=True):
        # NOTE: If drop_objects is True, then anything object column not in the index_col_name will be dropped
        assert os.path.isfile(filepath)
        if not isinstance(filepath, Path): filepath = Path(filepath)
        df = pd.read_csv(filepath, sep=sep)
        # Set index column
        df = df.set_index(self.index_col_name)
        if drop_objects: df = df.drop(columns=list(df.select_dtypes(include=['object']).columns))
        # NOTE: will always drop the first block of column name i.e. 'drop_keep1_keep2_...' --> 'keep1_keep2_...'
        df = df.rename(columns={x: "_".join(x.split("_")[1:]) for x in df.select_dtypes(include=['number']).columns})
        self.data_sources[filepath] = True

        # converting some columns to data indexes if necessary
        # for example, a df could have 6 columns: [n1, c1, n2, c2, n3, c3]
        # maybe [n1, c1, n2, c2] belong to one group while [n3, c3] belong to a second
        # this converts [n3, c3] --> [n1, c1] and merges DFs mith MultiIndex labels
        # TODO: verify this works for assymetric number of n, c for each merge.
        # TODO: either verify this works for >2 categories OR error check to limit to 2 categories
        if sample_number_cats is not None:
            dfs = {}
            for cat, numbers in sample_number_cats.items(): # loop through "cat": [n1, n2, ...]
                cols = [] # for identifying [n3, n4, ...]
                new_cols = [] # for replacing [n3, n4, ...] with [n1, n2, ...]
                for n in numbers: # loop through [n1, n2, ...]
                    for col in df.columns:
                        if str(n) in col: 
                            cols.append(col)
                            new_cols.append(col.replace(str(n), str(np.ceil(len(cols) / 2).astype(int))))
                dfs[cat] = df[cols]
                dfs[cat].columns = new_cols
            
            # DFs should now be split in 'dfs'
            df = pd.concat(dfs.values(), keys=dfs.keys(), names=[sample_number_cat_name, df.index.name])
            df = df.reorder_levels(df.index.names[::-1])

        # return DF
        return df
    
    def add_csv(self, filepath: PathLike, sep: str, 
                sample_number_cat_name: NoneType | str, sample_number_cats: NoneType | dict, 
                category_name, old_data_category, new_data_category, drop_objects: bool=True
        ):
        """
        Description: Will add a new df from a new filepath, and create a new column to distinguish
        between old and new data.

        NOTE: For all index columns, the same unique values must be present for old and new data
        for the proper merge to occur.

        Args:
        - category_name: name of the new column with categories separating old and new data
        - old_data_category: name of the category that old data will be filled with in the new column
        - new_data_category: name of the category that the new data will be filled with the new column
        """
        # TODO: add error checking for if indexes don't match between original and new df
        # TODO: decide whether to allow more than 2x2x2x... for each MultiIndex level. If 2x2 is maintained, then heatmaps always possible. If NxM categories allowed, then merging and heatmaps become trickier.
        add_df = self._load_csv(
            filepath, sep, 
            sample_number_cat_name=sample_number_cat_name, sample_number_cats=sample_number_cats,
            drop_objects=drop_objects
        )
        assert np.all(add_df.index == self.df.index), f"DF index must match, but {add_df.index} and {self.df.index} do not!"
        self.df = pd.concat(
            [self.df, add_df], keys=[old_data_category, new_data_category], 
            names=[category_name] + self.df.index.names
        )
        self.df = self.df.reorder_levels(
            self.df.index.names[1:] + [self.df.index.names[0]]
        )

    def _setup_comparison_table(self, index_level_name, filters: NoneType | dict):
        """
        Sets up comparison pivot table.

        Arg:
        - index_level_name: name of the MultiIndex level that will be used to index the pivot comparison table
        - filters: for other indexes, filter to limit the number of variables to 2.
        """
        # Validating inputs
        EXPECTED = 2
        if self.currently_comparing: raise ValueError(f"Already comparing {self.current_comparison}. Close current comparison before starting a new one!")
        this_index_level_n_unique = len(self.df.index.get_level_values(index_level_name).unique())
        assert this_index_level_n_unique == EXPECTED, f"Expected {EXPECTED} unique values in MultiIndex level {index_level_name}, got {this_index_level_n_unique}."
        if filters is not None:
            assert len(self.df.index.names) == len(filters) + 2, \
                f"DF MultiIndex has {len(self.df.index.names)} levels. Expected a filter of size {len(self.df.index.names) - 2}, but got {len(filters)}"

        # Creating MultiIndex based on filters
        for key in filters.keys(): assert key >= 0 and key < len(self.df.index.names), f"Filter key {key} not a valid level in DF MultiIndex."
        filter_slice = tuple(filters[idx] if idx in filters.keys() else slice(None) for idx in np.arange(len(self.df.index.names)))

        # Creating pivot table
        self.comparison_table = pd.pivot_table(
            data=self.df.loc[filter_slice,:],
            columns=index_level_name,
            values=self.df.columns,
            index=self.index_col_name,
        )
        self.currently_comparing = True
        self.current_comparison = [index_level_name, filters]

    def _comparison_cache_key(self, p_table, metadata_table, params: dict) -> str:
        """
        Content-addressed key: hashes the exact count matrix and metadata DESeq2 would see,
        plus every parameter that changes the result. If the data or settings change,
        the key changes, so stale results are never loaded.
        """
        h = hashlib.sha256()
        h.update(pd.util.hash_pandas_object(p_table, index=True).values.tobytes())
        h.update("\x1f".join(map(str, p_table.columns)).encode())  # gene order/names
        h.update(pd.util.hash_pandas_object(metadata_table, index=True).values.tobytes())
        h.update(json.dumps(params, sort_keys=True).encode())
        return h.hexdigest()[:20]

    def _cache_paths(self, key: str):
        return self.cache_dir / f"{key}.csv.gz", self.cache_dir / f"{key}.json"

    def _save_comparison_cache(self, key: str, params: dict):
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        res_path, meta_path = self._cache_paths(key)
        # Write to temp files, then rename, so an interrupted run can't leave a corrupt cache entry.
        tmp_res = res_path.with_name(res_path.name + ".tmp")
        self.comparison_results.to_csv(tmp_res, compression="gzip")
        os.replace(tmp_res, res_path)
    
        level, filters = self.current_comparison
        meta = {
            "level": level,
            "filters": {str(k): str(v) for k, v in (filters or {}).items()},
            "params": params,
            "n_genes": int(len(self.comparison_results)),
            "created": datetime.now().isoformat(timespec="seconds"),
        }
        tmp_meta = meta_path.with_name(meta_path.name + ".tmp")
        tmp_meta.write_text(json.dumps(meta, indent=2))
        os.replace(tmp_meta, meta_path)

    def list_cached_comparisons(self) -> pd.DataFrame:
        """Human-readable index of what's in the cache directory."""
        if self.cache_dir is None or not self.cache_dir.exists():
            return pd.DataFrame()
        rows = []
        for p in sorted(self.cache_dir.glob("*.json")):
            meta = json.loads(p.read_text())
            rows.append({"key": p.stem, "level": meta["level"], **{f"filter_{k}": v for k, v in meta["filters"].items()},
                        "n_genes": meta["n_genes"], "created": meta["created"]})
        return pd.DataFrame(rows)
        
    # ---------------------------------------------------------------------------
    # _calculate_comparison_p_values
    # ---------------------------------------------------------------------------
    def _calculate_comparison_p_values(self, n_cpus: int = 8, alpha: float = 0.05,
                                        use_cache: bool = True, force_recompute: bool = False):
        assert self.currently_comparing, "Must start comparing before calculating p-values!"
        level = self.current_comparison[0]
    
        count_cols = [c for c in self.comparison_table.columns.get_level_values(0).unique() if "count" in c]
        p_table = self.comparison_table.loc[:, count_cols].T  # rows: (replicate_col, condition), cols: genes
    
        # Traceable sample IDs instead of a bare RangeIndex, e.g. "1_count|Ctrl"
        sample_ids = ["|".join(map(str, t)) for t in p_table.index]
        metadata_table = pd.DataFrame({level: p_table.index.get_level_values(level)}, index=sample_ids)
        p_table.index = sample_ids
    
        # Round, don't truncate: astype(int) turns 0.9 into 0 and biases every gene downward.
        p_table = p_table.round().astype(int)
    
        # TODO: consider a fancier design like '~sample + treatment + sample:treatment'
        # Consider: https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
        # Consider: https://www.biostars.org/p/221122/
        # Consider: R-formula cheatsheet - https://www.econometrics.blog/post/the-r-formula-cheatsheet/
        design = f"~{level}"
        contrast = np.array([0, 1])
        params = {
            "design": design,
            "contrast": contrast.tolist(),
            "alpha": alpha,
            "refit_cooks": True,
            "cooks_filter": True,
            "independent_filter": True,
            "pydeseq2_version": pkg_version("pydeseq2"),
        }
    
        caching = use_cache and self.cache_dir is not None
        key = self._comparison_cache_key(p_table, metadata_table, params) if caching else None
    
        # --- Try the cache ---
        if caching and not force_recompute:
            res_path, _ = self._cache_paths(key)
            if res_path.exists():
                self.comparison_results = pd.read_csv(res_path, index_col=0)
                self.comparison_p = self.comparison_results["pvalue"].sort_values()
                self.dds, self.ds = None, None  # not cached; use force_recompute=True if you need them
                self.comparison_from_cache = True
                return
    
        # --- Compute ---
        inference = DefaultInference(n_cpus=n_cpus)
        self.dds = DeseqDataSet(
            counts=p_table,
            metadata=metadata_table,
            design=design,
            refit_cooks=True,
            inference=inference,
        )
        # Runs the canonical sequence, including fit_MAP_dispersions (dispersion shrinkage).
        self.dds.deseq2()
    
        # Reference level = alphabetically first condition unless you set ref_level,
        # so [0, 1] is (second condition) vs (first condition).
        self.ds = DeseqStats(
            self.dds,
            contrast=contrast,
            alpha=alpha,
            cooks_filter=True,
            independent_filter=True,
            inference=inference,
        )
        # summary() is where Cook's filtering, independent filtering, and BH adjustment happen.
        self.ds.summary()
        self.comparison_results = self.ds.results_df  # baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
        self.comparison_p = self.comparison_results["pvalue"].sort_values()
        self.comparison_from_cache = False
    
        if caching:
            self._save_comparison_cache(key, params)

    @staticmethod
    def _set_log_cpm_ticks(axis, lo: float, hi: float, pseudocount: float = 0.0,
                        label_style: str = "pow10", minor: bool = True):
        """
        Place ticks at raw CPM values (…, 0.1, 1, 10, 100, …, plus 2–9 × each decade as minor ticks),
        positioned at log10(CPM + pseudocount), and label them with the raw CPM value.
    
        - axis: ax.xaxis or ax.yaxis
        - lo, hi: axis limits in plotted (log10) units
        - label_style: "pow10" -> $10^{n}$, "plain" -> 0.1, 1, 10, 100
        """
        pc = pseudocount or 0.0
    
        # Raw-CPM decades to consider. With a pseudocount, values below pc are crushed
        # toward log10(pc), so start at pc's decade to avoid a pile of unreadable minor ticks.
        raw_hi = 10 ** hi - pc
        raw_lo = 10 ** lo - pc if pc == 0 else pc
        d_min = int(np.floor(np.log10(max(raw_lo, 1e-12))))
        d_max = int(np.ceil(np.log10(max(raw_hi, 1e-12))))
    
        def pos(v):
            return np.log10(v + pc)
    
        def in_range(p):
            return lo - 1e-9 <= p <= hi + 1e-9
    
        def fmt(d):
            if label_style == "plain":
                return f"{10.0 ** d:g}"
            return f"$10^{{{d}}}$"
    
        major_pos, major_lab = [], []
        if pc > 0 and in_range(pos(0.0)):
            major_pos.append(pos(0.0))
            major_lab.append("0")
        for d in range(d_min, d_max + 1):
            v = 10.0 ** d
            if pc > 0 and v < pc:
                continue
            p = pos(v)
            if in_range(p):
                major_pos.append(p)
                major_lab.append(fmt(d))
    
        minor_pos = []
        if minor:
            for d in range(d_min, d_max + 1):
                for m in range(2, 10):
                    v = m * 10.0 ** d
                    if pc > 0 and v < pc:
                        continue
                    p = pos(v)
                    if in_range(p):
                        minor_pos.append(p)
    
        axis.set_ticks(major_pos, labels=major_lab)
        axis.set_ticks(minor_pos, minor=True)
    
    
    # ---------------------------------------------------------------------------
    # 6) plot_comparison_density
    # ---------------------------------------------------------------------------
    def plot_comparison_density(
        self,
        ax=None,
        x_cond=None,
        y_cond=None,
        top_n: int | None = None,
        genes: list | dict | None=None,
        p_threshold: float | None = None,
        p_col: str = "padj",
        highlight_color: str | list="k",
        value_key: str = "cpm",
        pseudocount: float | None = None,
        min_log: float = -1.0,
        bins: int = 100,
        norm=None,
        cmap=None,
        colorbar: bool = True,
        log_ticks: bool = True,
        tick_label_style: str = "pow10",
        title: str | None = None,
        label_fontsize: int = 8,
        highlight_pointsize: int = 16,
        highlight_edgewidth: float = 1,
        max_labels: int = 50,
    ):
        """
        2D density of log10(mean CPM) for the current comparison, with genes highlighted by p-value.
    
        Highlight modes (mutually exclusive; none = no highlights):
        - top_n: the N genes with the lowest `p_col`
        - p_threshold: all genes with `p_col` <= threshold (labels capped at `max_labels`)
        - genes: an explicit list of gene names; their p-values are shown in the labels. Can also be
        a dictionary {"key1": gene_list1, "key2": gene_list2}
        - highlight_color: Either a single color (e.g. "k" for black) or a list. If a list, must
        be the same length as the genes dictionary.
    
        Args:
        - ax: existing Axes to draw into (for multi-panel figures). If None, a new figure is made.
        - x_cond, y_cond: condition labels for the axes. Default: alphabetical order, which
        matches the DESeq2 reference level used in _calculate_comparison_p_values.
        - p_col: "padj" (recommended) or "pvalue".
        - pseudocount: if None, genes with 0 CPM in either condition are dropped (log10(0) = -inf).
        Set e.g. 0.1 to keep them at the lower edge.
        - norm: pass a shared mcolors.LogNorm(vmin=1, vmax=...) to make colors comparable across panels.
        - log_ticks: label axes in raw CPM with log-style major/minor ticks (data stay in log10 units).
        - tick_label_style: "pow10" ($10^{n}$) or "plain" (0.1, 1, 10, ...).
    
        Returns: (ax, highlighted) where `highlighted` is a DataFrame of the labeled genes.
        """
        multiple_highlighted = False

        if not self.currently_comparing:
            raise ValueError("Start a comparison (_setup_comparison_table) before plotting.")
        if self.comparison_results is None:
            raise ValueError("Run _calculate_comparison_p_values() before plotting.")
        if sum(v is not None for v in (top_n, genes, p_threshold)) > 1:
            raise ValueError("Use only one of top_n, genes, p_threshold.")
        res = self.comparison_results
        if p_col not in res.columns:
            raise ValueError(f"p_col must be one of {list(res.columns)}, got {p_col!r}.")
    
        level, filters = self.current_comparison
    
        # --- Mean expression per condition across replicates ---
        tbl = self.comparison_table
        val_cols = [c for c in tbl.columns.get_level_values(0).unique() if value_key in c]
        if not val_cols:
            raise ValueError(f"No columns containing {value_key!r} in the comparison table.")
        expr = tbl.loc[:, val_cols].T.groupby(level=level).mean().T  # index: genes, cols: conditions
    
        conds = sorted(expr.columns)
        if len(conds) != 2:
            raise ValueError(f"Expected 2 conditions in level {level!r}, got {conds}.")
        x_cond = x_cond if x_cond is not None else conds[0]
        y_cond = y_cond if y_cond is not None else conds[1]
    
        # --- Log transform ---
        pc = 0.0 if pseudocount is None else pseudocount
        with np.errstate(divide="ignore"):
            plot_df = pd.DataFrame({
                "x": np.log10(expr[x_cond] + pc),
                "y": np.log10(expr[y_cond] + pc),
            })
        plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()
        plot_df = plot_df[(plot_df["x"] >= min_log) & (plot_df["y"] >= min_log)]
    
        # --- Select genes to highlight ---
        pvals = res[p_col]
        if top_n is not None:
            sel = pvals.dropna().nsmallest(top_n).index
        elif p_threshold is not None:
            sel = pvals[pvals <= p_threshold].sort_values().index
        elif genes is not None:
            if isinstance(genes, list): 
                sel = pd.Index(list(genes))
            if isinstance(genes, dict): 
                sel = {key: pd.Index(list(these_genes)) for key, these_genes in genes.items()}
                multiple_highlighted = True
        else:
            sel = pd.Index([])
    
        if multiple_highlighted:
            multi_sel = list(sel.values())[0]
            for s in list(sel.values())[1:]:
                multi_sel.union(s)
            not_in_results = multi_sel.difference(res.index)
        else:
            not_in_results = sel.difference(res.index)
        if len(not_in_results):
            warnings.warn(f"Not in DESeq2 results: {list(not_in_results)}")

        if multiple_highlighted:
            not_plotted = multi_sel.intersection(res.index).difference(plot_df.index)
        else:
            not_plotted = sel.intersection(res.index).difference(plot_df.index)
        if len(not_plotted):
            warnings.warn(
                f"{len(not_plotted)} highlighted gene(s) fall outside the plotted range "
                f"(zero CPM or below min_log): {list(not_plotted)}. Consider setting pseudocount."
            )
    
        if multiple_highlighted:
            df_list = []
            for key, s in sel.items():
                this_highlighted = plot_df.reindex(s).dropna()
                this_highlighted["key"] = key
                df_list.append(this_highlighted)
            highlighted = pd.concat(df_list, axis=0)
        else:
            highlighted = plot_df.reindex(sel).dropna()
        highlighted[p_col] = pvals.reindex(highlighted.index)
        highlighted = highlighted.sort_values(p_col, na_position="last")
        labeled = highlighted.head(max_labels)
        
        if len(highlighted) > max_labels:
            warnings.warn(f"{len(highlighted)} genes highlighted; labeling only the {max_labels} lowest {p_col}.")
    
        # --- Plot ---
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5.5), layout="constrained")
        else:
            fig = ax.figure
    
        if cmap is None:
            cmap = mpl.colormaps["rocket_r"].copy()  # copy so the global registry isn't modified
            cmap.set_bad(alpha=0)
        if norm is None:
            norm = mcolors.LogNorm()
    
        hi = float(np.ceil(max(plot_df["x"].max(), plot_df["y"].max()) * 10) / 10)
        lims = (min_log, hi)
        _, _, _, im = ax.hist2d(
            plot_df["x"], plot_df["y"],
            bins=bins, range=[lims, lims],
            cmap=cmap, norm=norm, cmin=1,
        )

        if multiple_highlighted:
            for color, key in zip(highlight_color, genes.keys()):
                ax.scatter(
                    highlighted[highlighted["key"]==key]["x"], 
                    highlighted[highlighted["key"]==key]["y"], 
                    color=color, s=highlight_pointsize, 
                    zorder=3, label=key, edgecolors="w", 
                    linewidth=highlight_edgewidth
                )
        else:
            ax.scatter(
                highlighted["x"], highlighted["y"], 
                color=highlight_color, s=highlight_pointsize, 
                zorder=3, edgecolors="w", 
                linewidth=highlight_edgewidth
            )

        if colorbar:
            fig.colorbar(im, ax=ax, label="genes per bin")

        ax.axline((0, 0), slope=1, color="0.5", lw=1, ls="--")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal")

        if multiple_highlighted: ax.legend()
    
        if log_ticks:
            self._set_log_cpm_ticks(ax.xaxis, *lims, pseudocount=pc, label_style=tick_label_style)
            self._set_log_cpm_ticks(ax.yaxis, *lims, pseudocount=pc, label_style=tick_label_style)
            unit = "mean CPM" + (f" (+{pc:g} pseudocount)" if pc else "")
        else:
            unit = "log10 mean CPM" + (f" (+{pc:g})" if pc else "")
    
        texts = []
        for gene, row in labeled.iterrows():
            p = row[p_col]
            if pd.isna(p): p_str = f"{gene} NA"
            elif p >= 0.05: p_str = f"{gene} NS"
            elif p == 0: p_str = f"{gene} p=0"
            else: p_str = f"{gene} p={p:.0e}"
            # p_str = "NA" if pd.isna(p) else f"{p:.0e}"
            texts.append(ax.text(
                row["x"], row["y"], p_str, fontsize=label_fontsize,
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
            ))
    
        if texts:
            fig.canvas.draw()  # finalize layout (aspect, colorbar, ticks) before adjusting
            adjust_text(
                texts,
                avoid_self=True,
                x=highlighted["x"].to_numpy(), y=highlighted["y"].to_numpy(),  # points to avoid
                ax=ax, expand=(1.2, 1.4), force_text=(0.5,0.5), force_static=(0.5,0.5),
                arrowprops=dict(arrowstyle="-", color="0.3", lw=1),
            )
    
        ax.set_xlabel(f"{x_cond} {unit}")
        ax.set_ylabel(f"{y_cond} {unit}")
        if title is None:
            context = " ".join(str(v) for v in filters.values()) if filters else ""
            if top_n is not None:
                desc = f"{top_n} lowest {p_col} highlighted"
            elif p_threshold is not None:
                desc = f"{p_col} \u2264 {p_threshold:g} highlighted"
            elif genes is not None:
                desc = "selected genes highlighted"
            else:
                desc = ""
            title = f"{context} CPM density: {y_cond} vs {x_cond}".strip() + (f"\n{desc}" if desc else "")
        ax.set_title(title)
    
        return ax, highlighted
    
    # ---------------------------------------------------------------------------
    # 3) close_comparison (needed to loop over comparisons for multi-panel figures)
    # ---------------------------------------------------------------------------
    def close_comparison(self):
        self.currently_comparing = False
        self.current_comparison = None
        self.comparison_table = None
        self.dds = None
        self.ds = None
        self.comparison_p = None
        self.comparison_results = None
        self.comparison_from_cache = False

    def head(self):
        return self.df.head()

    def info(self):
        print(f"{'Dataframe ID: ':<20}{self.base_id:>10}")
        print(f"{'Dataframe Rows: ':<20}{self.df.shape[0]:>10}")
        print(f"{'Dataframe Columns: ':<20}{self.df.shape[1]:>10}")
        print(f"{'Number of Samples: ':<20}{self.n_samples:>10}\n")
        print(f"{'Sample name: ':<20}{'Col Index':>10}")
        for idx in range(len(self.sample_names)):
            print(f"{self.sample_names[idx]:<20}{self.sample_col_indexes[idx]:>10}")

