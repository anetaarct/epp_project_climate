import math
import os
import random


os.makedirs("figures", exist_ok=True)
random.seed(42)

models = [
    ("Absolute climate", ["Temperature", "Precipitation", "DTR", "NDVI", "Absolute latitude", "Distance to coast", "Nest box"]),
    ("Short anomalies", ["Temperature anomaly", "Precipitation anomaly", "DTR", "NDVI", "Absolute latitude", "Distance to coast", "Nest box"]),
    ("Long anomalies", ["Temperature anomaly", "Precipitation anomaly", "DTR", "NDVI", "Absolute latitude", "Distance to coast", "Nest box"]),
    ("Daily weather variability", ["Daily temperature SD", "Daily precipitation SD", "Wet days > 1 mm", "NDVI", "Absolute latitude", "Distance to coast", "Nest box"]),
]

width = 980
left = 260
right = 60
top = 70
panel_gap = 52
row_h = 28
panel_title_h = 28
x_min, x_max = -1.2, 1.2
plot_w = width - left - right

height = top + sum(panel_title_h + row_h * len(p[1]) for p in models) + panel_gap * (len(models) - 1) + 70


def x_pos(x):
    return left + (x - x_min) / (x_max - x_min) * plot_w


svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
svg.append('<rect width="100%" height="100%" fill="#ffffff"/>')
svg.append('<text x="40" y="36" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1a2e40">Mock preview: environmental effects forest plot</text>')
svg.append('<text x="40" y="58" font-family="Arial, sans-serif" font-size="12" fill="#5c6670">Example values only. Layout matches the planned figure style.</text>')

y_cursor = top

for model_name, labels in models:
    svg.append(f'<text x="40" y="{y_cursor + 18}" font-family="Arial, sans-serif" font-size="15" font-weight="700" fill="#1a2e40">{model_name}</text>')
    panel_y0 = y_cursor + panel_title_h
    panel_h = row_h * len(labels)

    zero_x = x_pos(0)
    svg.append(f'<line x1="{zero_x:.1f}" y1="{panel_y0 - 6}" x2="{zero_x:.1f}" y2="{panel_y0 + panel_h - 10}" stroke="#555" stroke-width="1.2" stroke-dasharray="5,5"/>')

    for tick in [-1.0, -0.5, 0.0, 0.5, 1.0]:
        tx = x_pos(tick)
        svg.append(f'<line x1="{tx:.1f}" y1="{panel_y0 - 6}" x2="{tx:.1f}" y2="{panel_y0 + panel_h - 10}" stroke="#eeeeee" stroke-width="1"/>')
        if model_name == models[-1][0]:
            svg.append(f'<text x="{tx:.1f}" y="{height - 42}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" fill="#5c6670">{tick:g}</text>')

    for i, label in enumerate(reversed(labels)):
        y = panel_y0 + i * row_h + 11
        estimate = random.gauss(0, 0.28)
        half = random.uniform(0.22, 0.55)
        lower = max(x_min, estimate - half)
        upper = min(x_max, estimate + half)
        credible = lower > 0 or upper < 0
        color = "#1b9e77" if credible else "#666666"
        fill = "#1b9e77" if credible else "#ffffff"

        svg.append(f'<text x="235" y="{y + 4}" text-anchor="end" font-family="Arial, sans-serif" font-size="12" fill="#263238">{label}</text>')
        svg.append(f'<line x1="{x_pos(lower):.1f}" y1="{y}" x2="{x_pos(upper):.1f}" y2="{y}" stroke="{color}" stroke-width="2.2"/>')
        svg.append(f'<circle cx="{x_pos(estimate):.1f}" cy="{y}" r="5.2" fill="{fill}" stroke="#000000" stroke-width="0.8"/>')

    y_cursor += panel_title_h + panel_h + panel_gap

svg.append(f'<text x="{left + plot_w / 2:.1f}" y="{height - 15}" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#263238">Posterior effect on EPP probability (logit scale)</text>')

legend_y = 42
svg.append(f'<line x1="{width - 315}" y1="{legend_y}" x2="{width - 280}" y2="{legend_y}" stroke="#1b9e77" stroke-width="2.2"/>')
svg.append(f'<circle cx="{width - 297}" cy="{legend_y}" r="5.2" fill="#1b9e77" stroke="#000"/>')
svg.append(f'<text x="{width - 270}" y="{legend_y + 4}" font-family="Arial, sans-serif" font-size="11" fill="#263238">95% CrI excludes zero</text>')
svg.append(f'<line x1="{width - 315}" y1="{legend_y + 20}" x2="{width - 280}" y2="{legend_y + 20}" stroke="#666666" stroke-width="2.2"/>')
svg.append(f'<circle cx="{width - 297}" cy="{legend_y + 20}" r="5.2" fill="#ffffff" stroke="#000"/>')
svg.append(f'<text x="{width - 270}" y="{legend_y + 24}" font-family="Arial, sans-serif" font-size="11" fill="#263238">95% CrI overlaps zero</text>')

svg.append("</svg>")

with open("figures/forest_plot_mock_preview.svg", "w", encoding="utf-8") as f:
    f.write("\n".join(svg))

print("figures/forest_plot_mock_preview.svg")
