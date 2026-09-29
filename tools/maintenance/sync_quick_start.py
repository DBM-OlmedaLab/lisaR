#!/usr/bin/env python3
"""Synchronise the published Riaz ON snippets from the installed R script."""
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
START = '<!-- BEGIN RIAZ ON EXAMPLE -->'
END = '<!-- END RIAZ ON EXAMPLE -->'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    code = (ROOT / 'inst/examples/riaz-on-quick-start.R').read_text().rstrip()
    for relative, fence in [
        ('README.md', '```r'),
        ('vignettes/riaz-on-quick-start.Rmd', '```{r riaz-on, eval=FALSE}'),
    ]:
        path = ROOT / relative
        text = path.read_text()
        if text.count(START) != 1 or text.count(END) != 1:
            raise SystemExit(f'{relative}: expected one marked example')
        before, rest = text.split(START)
        _, after = rest.split(END)
        updated = before + START + '\n' + fence + '\n' + code + '\n```\n' + END + after
        if args.check and text != updated:
            raise SystemExit(f'{relative}: run tools/maintenance/sync_quick_start.py')
        if not args.check:
            path.write_text(updated)
    print('Riaz ON snippets match the canonical installed script.')


if __name__ == '__main__':
    main()
