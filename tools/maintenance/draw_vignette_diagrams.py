#!/usr/bin/env python3
"""Export documentation-only SVG and PNG from the same labelled diagram.

Requires Pillow only for the PNG fallback. Not a runtime/package dependency.
No scientific data are calculated or changed.
"""
from pathlib import Path
from html import escape
import math
import sys
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "vignettes" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
INK = "#183347"
TEAL = "#087E8B"
ONLY = set(sys.argv[1:])

class Diagram:
    def __init__(self, width, height, title, description):
        self.width, self.height = width, height
        self.image = Image.new("RGB", (width*2, height*2), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
                    f'<title id="title">{escape(title)}</title><desc id="desc">{escape(description)}</desc>',
                    '<rect width="100%" height="100%" fill="white"/>']
    def box(self, x, y, w, h, fill="#EEF5F8", stroke="#B9CBD5"):
        self.svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{fill}" stroke="{stroke}"/>')
        self.draw.rounded_rectangle((x*2,y*2,(x+w)*2,(y+h)*2), radius=20, fill=fill, outline=stroke, width=2)
    def text(self, x, y, value, size=20, bold=False, colour=INK):
        self.svg.append(f'<text x="{x}" y="{y}" font-family="DejaVu Sans, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{colour}">{escape(value)}</text>')
        self.draw.text((x*2,y*2), value, font=ImageFont.truetype(BOLD if bold else FONT,size*2), fill=colour, anchor="ls")
    def lines(self,x,y,values,size=18,step=27):
        for i,t in enumerate(values): self.text(x,y+i*step,t,size)
    def arrow(self,x1,y1,x2,y2,colour=TEAL):
        angle=math.atan2(y2-y1,x2-x1)
        pts=[(x2,y2),(x2-11*math.cos(angle-.45),y2-11*math.sin(angle-.45)),(x2-11*math.cos(angle+.45),y2-11*math.sin(angle+.45))]
        self.svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colour}" stroke-width="2.5"/>')
        self.svg.append('<polygon points="'+' '.join(f'{x},{y}' for x,y in pts)+f'" fill="{colour}"/>')
        self.draw.line((x1*2,y1*2,x2*2,y2*2), fill=colour,width=5)
        self.draw.polygon([(x*2,y*2) for x,y in pts],fill=colour)
    def save(self,name):
        if ONLY and name not in ONLY:
            return
        (OUT/f"{name}.svg").write_text("\n".join(self.svg+["</svg>\n"]))
        self.image.save(OUT/f"{name}.png",optimize=True)
    def connector(self, points):
        self.svg.append('<polyline points="'+' '.join(f'{x},{y}' for x,y in points[:-1])+f'" fill="none" stroke="{TEAL}" stroke-width="2.5"/>')
        self.draw.line([(x*2,y*2) for x,y in points[:-1]], fill=TEAL, width=5)
        self.arrow(*points[-2], *points[-1])

d=Diagram(1100,970,"From a prepared study to reusable results", "An upstream model supplies differential-expression results. lisaR runs enrichment, annotates gene sets with a pre-built dictionary and hierarchy, and saves the results and HTML report. Saved results branch into static browsing, selected figure generation and reuse, and portable export. These later tasks do not repeat DE, GSEA or ORA or change fixed assignments. Exported static reports contain no live generation controls.")
d.text(35,45,"From a prepared study to reusable results",28,True)
d.text(35,80,"Existing results enter the pipeline; semantic resources are prepared beforehand.",17)
d.box(25,220,218,165,"#F5F3EF")
d.text(42,251,"UPSTREAM",14,True)
d.text(42,283,"DE results",20,True)
d.lines(42,314,["Complete result table","Gene IDs and statistic","Effect and adjusted P"],16,25)
d.box(268,125,805,305,"#F7FBFC", "#C2DFE2")
d.text(287,155,"lisaR",19,True,TEAL)
d.box(287,190,220,66)
d.lines(302,216,["Gene-set memberships","TERM2GENE"],16,23)
d.box(287,310,220,95)
d.text(302,339,"Enrichment",19,True)
d.lines(302,367,["Ranked genes → GSEA","ORA: optional separate test"],14,23)
d.arrow(243,344,287,344)
d.arrow(394,256,394,310)
d.box(542,190,242,66)
d.lines(557,216,["LISA dictionary + hierarchy","Gene sets → categories"],15,23)
d.box(542,310,242,95)
d.text(557,338,"LISA annotation",19,True)
d.lines(557,366,["Category summaries","Gene-set and gene evidence"],15,23)
d.arrow(663,256,663,310)
d.arrow(507,358,542,358)
d.box(820,310,231,95)
d.text(836,338,"Saved results",19,True)
d.lines(836,366,["HTML report and tables","Inputs, settings and resources"],14,23)
d.arrow(784,358,820,358)
d.connector([(937,405),(937,477),(550,477),(550,530)])
d.text(35,488,"Once the run is saved",24,True)
d.text(35,518,"Choose the task you need.",17)
for x,title,lines in [
    (35,"Browse the report",["Open report_index.html","Inspect figures and tables","No R or Shiny needed"]),
    (395,"Generate selected figures",["Use Shiny or the R API","Reuse the saved results","Download data and scripts"]),
    (755,"Export and share",["Include available figures","Keep every subfolder","Open without R or Shiny"])]:
    d.box(x,555,310,172,"#EAF6F4")
    d.text(x+15,590,title,18,True)
    d.lines(x+15,629,lines,16,29)
d.connector([(550,530),(190,530),(190,555)])
d.arrow(550,530,550,555)
d.connector([(550,530),(910,530),(910,555)])
d.arrow(705,649,755,649)
d.text(35,772,"Browsing, figure generation and export do not repeat DE, GSEA or ORA.",18,True)
d.text(35,804,"Original scientific results and fixed LISA assignments stay unchanged.",18)
d.text(35,836,"Exported static copies have no live generation controls.",18)
d.text(35,878,"Optional profile contrast: A minus B, with each side's evidence retained.",16)
d.text(35,905,"Category summaries and profile differences are descriptive, not new tests.",16)
d.text(35,941,"No LLM is called when lisaR analyses a study.",18,True)
d.save("lisa-workflow")

d=Diagram(1100,470,"The resources used by lisaR", "Memberships link genes to gene sets, the dictionary links gene sets to categories with tier, and the category map places categories within supercategories. A native KEGG map is a separate optional graphical resource.")
d.text(35,45,"Resources have different roles",28,True)
columns=[(35,"Gene-set memberships",["Gene set ↔ gene","Defines the sets used by GSEA","and optional ORA."],"Obtained separately"),(395,"LISA dictionary",["Gene set ↔ category","Stores accepted assignments","and their tier."],"Core and expanded included"),(755,"Category hierarchy",["Category → supercategory","Organises related categories","for summaries and navigation."],"Included")]
for x,title,lines,foot in columns:
    d.box(x,110,310,235)
    d.text(x+16,145,title,19,True)
    d.lines(x+16,193,lines,16,30)
    d.text(x+16,308,foot,15,True,TEAL)
d.arrow(348,214,390,214); d.arrow(708,214,750,214)
d.text(35,398,"Native KEGG maps are optional graphics, separate from KEGG gene-set memberships.",18)
d.text(35,432,"Pre-built annotations are reused locally; they are not inferred from the study.",18)
d.save("lisa-resources")

d=Diagram(1060,560,"Riaz selected FULL design", "Two differential-expression comparisons, responders versus progressive disease before and during nivolumab, feed two LISA profiles. The descriptive contrast is during minus before. The standard example retains five analyses and two contrasts.")
d.text(35,45,"Riaz: the selected FULL comparison",28,True)
for x,letter,visit,ident in [(35,"B","Before treatment","responders_vs_pd_pre"),(570,"A","During treatment","responders_vs_pd_on")]:
    d.box(x,103,450,190)
    d.text(x+20,137,f"{letter} · {visit}",22,True)
    d.lines(x+20,181,["Responders versus progressive disease","DESeq2 Wald statistic → LISA profile"],18,32)
    d.text(x+20,264,ident,17,False,TEAL)
d.arrow(260,293,400,365); d.arrow(795,293,660,365)
d.box(250,365,560,106,"#EAF6F4")
d.text(278,404,"Profile contrast: A minus B",23,True)
d.text(278,440,"During minus before; a descriptive comparison",18)
d.text(35,518,"STANDARD: 5 DE analyses + 2 contrasts. Selected FULL: 2 DE analyses + 1 contrast.",18)
d.save("riaz-design")

d=Diagram(1100,460,"CPTAC ccRCC prepared input", "The prepared example has eighty matched tumour and adjacent tissue pairs, or one hundred and sixty samples. Upstream paired modelling and protein-to-gene mapping produce 6482 unique gene rows for one differential-abundance comparison. lisaR analyses four semantic collections. The saved audit does not give the total originally measured protein-feature count.")
d.text(35,45,"CPTAC ccRCC: from paired samples to a gene-level input",26,True)
for x, title, lines in [(35,"80 matched pairs",["Tumour + adjacent tissue","160 samples","Paired differential abundance"]), (400,"6,482 gene rows",["Mapped and filtered upstream","Unique gene symbols","Signed limma t for ranking"]), (765,"One LISA analysis",["Tumour minus adjacent tissue","Four semantic collections","Category and gene evidence"])]:
    d.box(x,115,300,210)
    d.text(x+15,154,title,20,True)
    d.lines(x+15,203,lines,16,32)
d.arrow(339,219,395,219)
d.arrow(704,219,760,219)
d.text(35,384,"6,482 describes the prepared gene table, not all measured protein features.",18)
d.text(35,421,"The saved audit has no pre-mapping denominator for a mapping percentage.",18)
d.save("cptac-coverage")
