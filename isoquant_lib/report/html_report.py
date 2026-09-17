############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Renders a RunSummary as a single self-contained HTML page.

No templating engine and no new dependency: the page is plain HTML with inline CSS,
and the figures are matplotlib SVGs embedded directly in the document, so the file
can be copied or mailed around on its own. If matplotlib is unavailable the tables
are still written - only the figures are dropped.
"""

import contextlib
import html
import io
import logging
from typing import List, Optional, Tuple

logger = logging.getLogger('IsoQuant')


def _pyplot():
    """Import pyplot for the report figures, or return None when it is unusable.

    Figures are optional: an environment with a broken or missing matplotlib still
    gets the tables. The import is done with stderr captured because a matplotlib
    built against a different numpy prints a traceback of its own on import.
    """
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        return plt
    except Exception as e:
        logger.debug("Report figures disabled, matplotlib is not usable: %s" % e)
        return None

PAGE_CSS = """
:root { color-scheme: light dark; }
body { font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
       margin: 0 auto; max-width: 62rem; padding: 2rem 1.25rem 4rem; line-height: 1.5;
       color: #1c1f24; background: #ffffff; }
h1 { font-size: 1.6rem; margin: 0 0 .25rem; }
h2 { font-size: 1.15rem; margin: 2.2rem 0 .75rem; padding-bottom: .3rem;
     border-bottom: 1px solid #e3e6ea; }
.subtitle { color: #6b7280; margin: 0 0 1.5rem; font-size: .9rem; word-break: break-all; }
.tiles { display: flex; flex-wrap: wrap; gap: .75rem; margin: 1rem 0 0; }
.tile { flex: 1 1 10rem; border: 1px solid #e3e6ea; border-radius: .5rem; padding: .75rem .9rem; }
.tile .value { font-size: 1.45rem; font-weight: 600; }
.tile .label { color: #6b7280; font-size: .8rem; text-transform: uppercase;
               letter-spacing: .03em; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; }
th, td { text-align: left; padding: .35rem .6rem; border-bottom: 1px solid #eef0f3; }
th { color: #6b7280; font-weight: 600; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.table-wrap { overflow-x: auto; }
figure { margin: 1rem 0 0; }
figure svg { max-width: 100%; height: auto; }
footer { margin-top: 3rem; color: #6b7280; font-size: .8rem; }
@media (prefers-color-scheme: dark) {
  body { color: #e6e8eb; background: #16181c; }
  h2 { border-color: #2a2e35; }
  .tile { border-color: #2a2e35; }
  th, td { border-color: #24272d; }
}
"""


def _format_number(value) -> str:
    if value is None:
        return "&mdash;"
    if isinstance(value, float):
        if value.is_integer():
            return "{:,}".format(int(value))
        return "{:,.2f}".format(value)
    if isinstance(value, int):
        return "{:,}".format(value)
    return html.escape(str(value))


def _format_percent(rate) -> str:
    if rate is None:
        return "&mdash;"
    return "{:.1f}%".format(rate * 100.0)


def _tile(label: str, value: str) -> str:
    return ('<div class="tile"><div class="value">%s</div>'
            '<div class="label">%s</div></div>' % (value, html.escape(label)))


def _table(rows: List[Tuple[str, str]], headers: Tuple[str, str] = ("", "Reads")) -> str:
    if not rows:
        return ""
    body = "".join('<tr><td>%s</td><td class="num">%s</td></tr>' % (html.escape(str(name)), value)
                   for name, value in rows)
    return ('<div class="table-wrap"><table><thead><tr><th>%s</th><th class="num">%s</th></tr>'
            '</thead><tbody>%s</tbody></table></div>'
            % (html.escape(headers[0]), html.escape(headers[1]), body))


def _bar_chart_svg(title: str, labels: List[str], values: List[float]) -> str:
    """Horizontal bar chart as an inline SVG, or "" when matplotlib is unusable."""
    if not labels:
        return ""
    plt = _pyplot()
    if plt is None:
        return ""
    try:
        height = max(1.6, 0.42 * len(labels) + 0.9)
        figure, axes = plt.subplots(figsize=(7.5, height))
        positions = range(len(labels))
        axes.barh(list(positions), values, color="#4c78a8")
        axes.set_yticks(list(positions))
        axes.set_yticklabels(labels, fontsize=9)
        axes.invert_yaxis()
        axes.set_title(title, fontsize=10, loc="left")
        axes.spines["top"].set_visible(False)
        axes.spines["right"].set_visible(False)
        axes.tick_params(axis="x", labelsize=8)
        figure.tight_layout()
        buffer = io.StringIO()
        figure.savefig(buffer, format="svg")
        plt.close(figure)
    except Exception as e:
        logger.debug("Cannot render report figure '%s': %s" % (title, e))
        return ""

    svg = buffer.getvalue()
    # Drop the XML prolog and DOCTYPE so the SVG can be inlined in the page body.
    start = svg.find("<svg")
    return "<figure>%s</figure>" % svg[start:] if start >= 0 else ""


def _kpi_tiles(summary) -> str:
    tiles = []
    if summary.total_reads is not None:
        tiles.append(_tile("Input reads", _format_number(summary.total_reads)))
        tiles.append(_tile("Mapping rate", _format_percent(summary.mapping_rate)))
    rollup = summary.assignment_rollup()
    if rollup:
        tiles.append(_tile("Uniquely assigned",
                           _format_percent(_safe_rate(rollup.get("unique"), rollup.get("total")))))
    if summary.barcodes:
        tiles.append(_tile("Valid barcodes", _format_percent(summary.barcode_rate)))
    for strategy, per_feature in (summary.groups or {}).items():
        gene_stats = summary.group_rollup(strategy).get("gene") or {}
        if gene_stats.get("groups") is not None:
            tiles.append(_tile("Barcodes/spots with counts", _format_number(gene_stats["groups"])))
        if gene_stats.get("share_of_reads") is not None:
            tiles.append(_tile("Reads in barcodes (gene level)",
                               _format_percent(gene_stats["share_of_reads"])))
        break
    if not tiles:
        return ""
    return '<div class="tiles">%s</div>' % "".join(tiles)


def _safe_rate(numerator, denominator) -> Optional[float]:
    if not denominator or numerator is None:
        return None
    return numerator / denominator


def _alignment_section(summary) -> str:
    if not summary.alignment:
        return ""
    rows = [(name.replace("_", " ").capitalize(), _format_number(count))
            for name, count in summary.alignment.items()]
    rows.append(("Input reads (primary + unaligned)", _format_number(summary.total_reads)))
    rows.append(("Mapping rate", _format_percent(summary.mapping_rate)))
    return "<h2>Alignment</h2>" + _table(rows, ("Alignment records", "Count"))


def _assignment_section(summary) -> str:
    if not summary.assignment:
        return ""
    rollup = summary.assignment_rollup()
    rows = [(name.replace("_", " "), _format_number(count))
            for name, count in summary.assignment.items()]
    total = rollup.get("total")
    for label in ("unique", "ambiguous", "inconsistent", "unassigned"):
        rows.append(("%s (share)" % label.capitalize(),
                     _format_percent(_safe_rate(rollup.get(label), total))))
    if summary.total_assignments is not None:
        rows.append(("PolyA tail detected",
                     "%s (%s)" % (_format_number(summary.polya_reads),
                                  _format_percent(_safe_rate(summary.polya_reads,
                                                             summary.total_assignments)))))
    figure = _bar_chart_svg("Reads per assignment type",
                            [name.replace("_", " ") for name in summary.assignment],
                            list(summary.assignment.values()))
    return "<h2>Read assignment</h2>" + figure + _table(rows, ("Assignment type", "Reads"))


def _quantification_section(summary) -> str:
    if not summary.quantification:
        return ""
    rows = []
    for feature, stats in summary.quantification.items():
        name = feature.capitalize()
        counted = stats.get("counted")
        rows.append(("%s: counted reads" % name, _format_number(counted)))
        rows.append(("%s: share of input reads" % name,
                     _format_percent(_safe_rate(counted, summary.total_reads))))
        rows.append(("%s: ambiguous" % name, _format_number(stats.get("ambiguous"))))
        rows.append(("%s: no feature" % name, _format_number(stats.get("no_feature"))))
    return "<h2>Quantification</h2>" + _table(rows, ("Feature", "Reads"))


def _models_section(summary) -> str:
    if not summary.transcript_models:
        return ""
    rows = [(name.replace("_", " "), _format_number(count))
            for name, count in summary.transcript_models.items()]
    rows.append(("Total", _format_number(sum(summary.transcript_models.values()))))
    return "<h2>Discovered transcript models</h2>" + _table(rows, ("Model type", "Transcripts"))


def _barcode_section(summary) -> str:
    if not summary.barcodes and not summary.cell_barcodes and not summary.umi_filtering:
        return ""
    rows = [(name, _format_number(count)) for name, count in summary.barcodes.items()]
    if summary.barcodes:
        rows.append(("Valid barcode rate", _format_percent(summary.barcode_rate)))
    for name, count in summary.cell_barcodes.items():
        rows.append((name, _format_number(count)))
    return "<h2>Barcode calling</h2>" + _table(rows, ("Stage", "Reads"))


def _umi_section(summary) -> str:
    if not summary.umi_filtering:
        return ""
    rows = [(name, _format_number(count)) for name, count in summary.umi_filtering.items()]
    saved = summary.umi_filtering.get("Total reads saved")
    processed = summary.umi_filtering.get("Total assignments processed")
    if saved is not None and processed:
        rows.append(("Duplication rate", _format_percent((processed - saved) / processed)))
    title = "UMI deduplication"
    if summary.umi_edit_distance is not None:
        title += " (edit distance %d)" % summary.umi_edit_distance
    return "<h2>%s</h2>" % html.escape(title) + _table(rows, ("Stage", "Reads"))


def _rank_plot_svg(title: str, ranked_reads) -> str:
    """Barcode rank curve (depth vs rank, log-log), as an inline SVG."""
    if not ranked_reads:
        return ""
    plt = _pyplot()
    if plt is None:
        return ""
    try:
        figure, axes = plt.subplots(figsize=(7.5, 3.2))
        axes.plot(range(1, len(ranked_reads) + 1), ranked_reads, color="#4c78a8")
        axes.set_xscale("log")
        axes.set_yscale("log")
        axes.set_xlabel("Barcode rank", fontsize=9)
        axes.set_ylabel("Reads", fontsize=9)
        axes.set_title(title, fontsize=10, loc="left")
        axes.spines["top"].set_visible(False)
        axes.spines["right"].set_visible(False)
        axes.tick_params(labelsize=8)
        figure.tight_layout()
        buffer = io.StringIO()
        figure.savefig(buffer, format="svg")
        plt.close(figure)
    except Exception as e:
        logger.debug("Cannot render report figure '%s': %s" % (title, e))
        return ""
    svg = buffer.getvalue()
    start = svg.find("<svg")
    return "<figure>%s</figure>" % svg[start:] if start >= 0 else ""


def _groups_section(summary) -> str:
    if not summary.groups:
        return ""
    blocks = []
    for strategy in summary.groups:
        rollup = summary.group_rollup(strategy)
        gene_stats = (summary.groups.get(strategy) or {}).get("gene") or {}
        figure = _rank_plot_svg("Reads per group, ranked (%s)" % strategy,
                                gene_stats.get("ranked_reads"))
        rows = []
        for feature, stats in rollup.items():
            name = feature.capitalize()
            rows.append(("%s: groups with counts" % name, _format_number(stats.get("groups"))))
            rows.append(("%s: reads in groups" % name, _format_number(stats.get("reads"))))
            rows.append(("%s: share of input reads" % name,
                         _format_percent(stats.get("share_of_reads"))))
            rows.append(("%s: median reads per group" % name,
                         _format_number(stats.get("median_reads_per_group"))))
            rows.append(("%s: median features per group" % name,
                         _format_number(stats.get("median_features_per_group"))))
        blocks.append("<h2>Grouped counts: %s</h2>" % html.escape(strategy)
                      + figure + _table(rows, ("Metric", "Value")))
    return "".join(blocks)


def render_html(summary, file_name: str) -> None:
    """Write the summary as a standalone HTML page."""
    sections = [
        _kpi_tiles(summary),
        _alignment_section(summary),
        _assignment_section(summary),
        _quantification_section(summary),
        _models_section(summary),
        _barcode_section(summary),
        _umi_section(summary),
        _groups_section(summary),
    ]
    subtitle = "IsoQuant %s" % summary.isoquant_version if summary.isoquant_version else "IsoQuant"
    if summary.mode:
        subtitle += " &middot; mode: %s" % html.escape(str(summary.mode))
    page = (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<title>IsoQuant summary: %s</title>\n<style>%s</style>\n</head>\n<body>\n"
        "<h1>%s</h1>\n<p class=\"subtitle\">%s</p>\n%s\n"
        "<footer>%s</footer>\n</body>\n</html>\n"
        % (html.escape(summary.sample_name), PAGE_CSS, html.escape(summary.sample_name),
           subtitle, "\n".join(s for s in sections if s),
           html.escape(summary.command_line))
    )
    with open(file_name, "w") as f:
        f.write(page)
