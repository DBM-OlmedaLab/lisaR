gene_evidence_fixture <- function() {
  memberships <- data.frame(gs_name = c("SET_A", "SET_A", "SET_B", "SET_B", "SET_C", "SET_C", "SET_D"),
    gene_symbol = c("FAV", "OTHER", "FAV", "OTHER", "FAV", "MISSING_DE", "FAV"))
  ledger <- data.frame(pathway = c("SET_A", "SET_B", "SET_C", "SET_D"),
    NES = c(2, -1.4, .8, NA_real_), padj = c(.01, .02, .8, NA_real_),
    pval = c(.002, .006, .3, NA_real_), leadingEdge = c("OTHER", "FAV", "MISSING_DE", NA),
    eligible_for_gsea = c(TRUE, TRUE, TRUE, FALSE), gsea_eligibility_status = c(rep("eligible", 3), "below_min_size"),
    classification_status = c("classified", "unclassified", "classified", "classified"), analysis_collection = "GOBP-C2")
  annotation <- ledger; annotation$category_id <- c("CAT", NA, "CAT", "CAT")
  de <- data.frame(symbol = c("FAV", "OTHER"), gene_id = c("00001", "00002"),
    log2FoldChange = c(-.8, 2), rank_value = c(-3, 5), pvalue = c(.1, .0001), padj = c(.25, .001))
  list(scopes = list(list(analysis_id = "A", collection = "GOBP-C2", positive_contrast = "treated minus control",
    universe_ledger = ledger, gsea_table = annotation, de_table = de)), memberships = memberships)
}

test_that("global evidence uses full membership, not leading edges, and preserves unavailable states", {
  f <- gene_evidence_fixture()
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  q <- lisaR:::query_lisa_gene_evidence(e, "00001")
  expect_identical(q$status, "matched"); expect_identical(q$symbol, "FAV")
  expect_setequal(q$sets$pathway, c("SET_A", "SET_B"))
  expect_identical(q$sets$leading_edge[q$sets$pathway == "SET_A"], "no")
  expect_identical(q$sets$leading_edge[q$sets$pathway == "SET_B"], "yes")
  expect_identical(q$sets$classification_status[q$sets$pathway == "SET_B"], "unclassified")
  expect_equal(q$de$de_fdr, .25); expect_true(all(is.na(q$de$statistic)))
  all <- lisaR:::query_lisa_gene_evidence(e, "FAV", include_all = TRUE)
  expect_setequal(all$sets$pathway, c("SET_A", "SET_B", "SET_C", "SET_D"))
  expect_identical(all$sets$leading_edge[all$sets$pathway == "SET_D"], "unavailable")
  expect_true(is.na(all$sets$NES[all$sets$pathway == "SET_D"]))
  missing <- lisaR:::query_lisa_gene_evidence(e, "MISSING_DE", include_all = TRUE)
  expect_identical(missing$status, "matched"); expect_equal(nrow(missing$de), 0)
  expect_false(any(grepl("priorit|score", names(e$de))))
})

test_that("gene queries preserve numeric-looking IDs and never resolve ambiguous names silently", {
  f <- gene_evidence_fixture()
  f$scopes[[1]]$de_table$gene_id <- c("00001", "FAV")
  duplicate <- f$scopes[[1]]$de_table[1, ]; duplicate$gene_id <- "00003"
  f$scopes[[1]]$de_table <- rbind(f$scopes[[1]]$de_table, duplicate)
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  q <- lisaR:::query_lisa_gene_evidence(e, "FAV")
  expect_identical(q$status, "ambiguous"); expect_setequal(q$candidates, c("FAV", "OTHER"))
  exact <- lisaR:::query_lisa_gene_evidence(e, "00001")
  expect_identical(exact$status, "matched"); expect_equal(nrow(exact$de), 2)
  expect_identical(lisaR:::query_lisa_gene_evidence(e, "fav")$status, "not_found")
  path <- tempfile(fileext = ".tsv")
  utils::write.table(f$scopes[[1]]$de_table, path, sep = "\t", quote = FALSE, row.names = FALSE)
  f$scopes[[1]]$de_table <- path
  from_file <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  expect_true("00001" %in% from_file$de$gene_id)
})

test_that("shared memberships and DE are not copied for each category or collection", {
  f <- gene_evidence_fixture(); second <- f$scopes[[1]]; second$collection <- "PATHWAYS"
  second$universe_ledger$analysis_collection <- "PATHWAYS"; second$gsea_table$analysis_collection <- "PATHWAYS"
  f$scopes[[2]] <- second
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  expect_equal(nrow(e$de), 2); expect_equal(nrow(e$memberships), 7)
  expect_equal(nrow(e$sets), 8); expect_equal(nrow(e$scopes), 2)
  f$scopes[[2]]$de_table$log2FoldChange[[1]] <- 0
  expect_error(lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core"), "DE evidence differs")
})

test_that("global gene evidence rejects wrong scope, corrupted membership and inconsistent results", {
  f <- gene_evidence_fixture(); x <- f
  x$scopes[[1]]$collection <- "GOMF"
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "collection scope")
  x <- f; x$scopes[[1]]$universe_ledger$analysis_id <- "OTHER_ANALYSIS"
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "analysis scope")
  x <- f; x$scopes[[1]]$gsea_table$tier <- "expanded"
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "tier scope")
  x <- f; x$scopes[[1]]$analysis_type <- "contrast"
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "individual analyses only")
  x <- f; x$memberships <- x$memberships[-2, ]
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "Leading-edge gene is not a member")
  x <- f; x$scopes[[1]]$gsea_table$NES[[1]] <- 10
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "ledger and annotation disagree")
  x <- f; row <- x$scopes[[1]]$gsea_table[1, ]; row$NES <- 10
  x$scopes[[1]]$gsea_table <- rbind(x$scopes[[1]]$gsea_table, row)
  expect_error(lisaR:::build_lisa_gene_evidence(x$scopes, x$memberships, "core"), "Conflicting GSEA")
})

test_that("full source memberships follow existing GSEA uppercase convention with spelling retained", {
  f <- gene_evidence_fixture(); f$memberships$gene_symbol <- tolower(f$memberships$gene_symbol)
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  q <- lisaR:::query_lisa_gene_evidence(e, "FAV", include_all = TRUE)
  expect_equal(nrow(q$sets), 4)
  expect_true("fav" %in% e$source_symbols$source_symbol)
  expect_true("FAV" %in% e$memberships$symbol)
  expect_false("fav" %in% e$memberships$symbol)
})

test_that("recorded DE provenance distinguishes row indexes and actual test-statistic rankings", {
  f <- gene_evidence_fixture()
  f$scopes[[1]]$column_mapping <- data.frame(standard_column = c("gene_id", "rank_value"), source_column = c("rowname", "stat"))
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  expect_identical(lisaR:::query_lisa_gene_evidence(e, "00001")$status, "not_found")
  q <- lisaR:::query_lisa_gene_evidence(e, "FAV")
  expect_equal(q$de$statistic, -3); expect_identical(q$de$gene_id_source, "rowname")
  f$scopes[[1]]$column_mapping$source_column[[2L]] <- "log2FoldChange"
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  expect_true(all(is.na(e$de$statistic)))
})

test_that("offline gene payload stays shared, exact and safe for data-bearing HTML", {
  f <- gene_evidence_fixture(); f$scopes[[1]]$positive_contrast <- "A </script> B"
  e <- lisaR:::build_lisa_gene_evidence(f$scopes, f$memberships, "core")
  payload <- lisaR:::lisa_gene_evidence_payload(e)
  expect_equal(length(payload$member_sets), length(unique(e$identifiers$symbol)))
  fav <- match("FAV", payload$symbols)
  expect_setequal(payload$pathways[payload$member_sets[[fav]] + 1L], c("SET_A", "SET_B", "SET_C", "SET_D"))
  root <- tempfile("gene-evidence-")
  out <- lisaR:::render_lisa_gene_evidence(e, root)
  expect_identical(out$status, "completed")
  html <- paste(readLines(out$html, warn = FALSE), collapse = "\n")
  expect_false(grepl("A </script> B", html, fixed = TRUE))
  expect_true(grepl("\\u003c", html, fixed = TRUE))
  lines <- readLines(out$html, warn = FALSE)
  payload_line <- lines[startsWith(lines, '<script type="application/json" id="gene-evidence-data">')]
  payload_line <- substring(payload_line, nchar('<script type="application/json" id="gene-evidence-data">') + 1L)
  payload_line <- substr(payload_line, 1L, nchar(payload_line) - nchar('</script>'))
  parsed <- jsonlite::fromJSON(payload_line)
  expect_identical(parsed$analyses$positive_contrast, "A </script> B")
  expect_true(file.exists(file.path(root, "assets", "gene-evidence.js")))
  expect_true(file.exists(file.path(root, "tables", "memberships.tsv")))
  expect_false(file.exists(file.path(root, "tables", "PATHWAYS_unclassified.tsv")))
  expect_match(html, 'id="scope-error" role="alert" hidden', fixed = TRUE)
  expect_match(html, 'id="reset-scope" href="index.html" hidden', fixed = TRUE)
  expect_match(html, 'id="contrast-return" hidden', fixed = TRUE)
})

test_that("gene URL scope and fixed local contrast return preserve exact identities", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; browser scope regression needs its URL API")
  script <- system.file("gene-evidence", "viewer.js", package = "lisaR")
  skip_if(!nzchar(script))
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "const source=fs.readFileSync(process.argv[2],'utf8');",
    "const start=source.indexOf('  function resolveGeneScope(query) {'),end=source.indexOf('  // End URL scope helpers.',start);",
    "assert.ok(start>=0&&end>start);",
    "const data={metadata:{tier:'core'},analyses:[{analysis_id:'A + 001'},{analysis_id:'B & 002'},{analysis_id:'C'}],scopes:[{analysis_id:'A + 001',collection:'GOBP-C2'},{analysis_id:'B & 002',collection:'GOBP-C2'},{analysis_id:'C',collection:'GOCC'}]};",
    "const api=vm.runInNewContext(source.slice(start,end)+';({resolveGeneScope,resolveContrastOrigin,contrastReturnHref})',{data,URLSearchParams});",
    "const initial=new URLSearchParams({analysis_id:'A + 001',collection:'GOBP-C2',tier:'core',gene:'001+A&B'});",
    "const single=api.resolveGeneScope(initial);assert.equal(single.error,'');assert.equal(single.analysis,'A + 001');assert.equal(single.collection,'GOBP-C2');assert.equal(api.resolveContrastOrigin(initial).origin,null);",
    "const all=api.resolveGeneScope(new URLSearchParams());assert.equal(all.error,'');assert.equal(all.analysis,'');assert.equal(all.collection,'');",
    "for(const [key,value]of[['tier','expanded'],['analysis_id','unknown'],['collection','unknown'],['collection','GOCC']]){const q=new URLSearchParams(initial);q.set(key,value);assert.ok(api.resolveGeneScope(q).error);}",
    "for(const side of ['a','b']){",
    "  const q=new URLSearchParams(initial);q.set('analysis_id',side==='a'?'A + 001':'B & 002');",
    "  Object.entries({return_contrast_scope:'C01_custom output',return_contrast_id:'C01',return_category_id:'0007',return_analysis_a:'A + 001',return_analysis_b:'B & 002',return_side:side}).forEach(([k,v])=>q.set(k,v));",
    "  const result=api.resolveContrastOrigin(q);assert.equal(result.error,'');assert.equal(result.origin.side,side);",
    "  const href=api.contrastReturnHref(result.origin,'002+NEW&GENE');assert.ok(href.startsWith('../contrast_evidence/'));",
    "  const url=new URL(href,'file:///D:/evaluation/report_pages/gene_evidence/index.html');",
    "  assert.equal(url.protocol,'file:');assert.equal(url.pathname,'/D:/evaluation/report_pages/contrast_evidence/C01_custom%20output/GOBP-C2/index.html');",
    "  assert.deepEqual(Object.fromEntries(new URLSearchParams(url.hash.slice(1))),{contrast_id:'C01',analysis_a:'A + 001',analysis_b:'B & 002',collection:'GOBP-C2',tier:'core',category_id:'0007',gene:'002+NEW&GENE'});",
    "  for(const [key,value]of[['return_contrast_scope','..'],['return_contrast_scope','https://outside.invalid/'],['return_contrast_scope','other/path'],['return_contrast_scope',''],['return_side','neither'],['return_analysis_b','unknown'],['return_analysis_b','A + 001'],['return_category_id',''],['tier','expanded'],['collection','GOCC']]){const bad=new URLSearchParams(q);bad.set(key,value);assert.ok(api.resolveContrastOrigin(bad).error);assert.equal(api.resolveContrastOrigin(bad).origin,null);}",
    "  const wrongSide=new URLSearchParams(q);wrongSide.set('analysis_id',side==='a'?'B & 002':'A + 001');assert.ok(api.resolveContrastOrigin(wrongSide).error);",
    "}",
    "process.stdout.write('PASS exact scope, recoverable rejection, and local A/B contrast return\\n');"
  ), probe, useBytes = TRUE)
  result <- suppressWarnings(system2(node, c(shQuote(probe), shQuote(script)), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS exact scope, recoverable rejection, and local A/B contrast return", fixed = TRUE)
})


test_that("gene filter changes survive reload without changing contrast return scope", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; filter URL regression needs its URL API")
  script <- system.file("gene-evidence", "viewer.js", package = "lisaR")
  skip_if(!nzchar(script))
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "const source=fs.readFileSync(process.argv[2],'utf8');",
    "const A='A + 001',B='B & 002',C='C';",
    "const data={metadata:{tier:'core',gsea_padj_cutoff:0.05,de_padj_cutoff:0.05},",
    "  symbols:['001+A&B'],pathways:['SET'],member_sets:[[0]],leading_edges:[[0],[],[0]],",
    "  analyses:[A,B,C].map(analysis_id=>({analysis_id,positive_contrast:'treated minus control'})),",
    "  scopes:[{analysis_id:A,collection:'GOBP-C2'},{analysis_id:B,collection:'GOBP-C2'},{analysis_id:C,collection:'GOCC'}],",
    "  de_columns:['symbol','gene_id','gene_id_source','analysis_id','log2FC','rank_value','statistic','pvalue','de_fdr'],",
    "  de:[[0],['00001'],['gene_id'],[A],[1.5],[2],[2],[0.01],[0.02]],",
    "  set_columns:['pathway','analysis_id','collection','NES','gsea_fdr','significant','leading_edge_state','evaluable','category_ids','eligibility_status'],",
    "  sets:[[0,0,0],[A,B,C],['GOBP-C2','GOBP-C2','GOCC'],[2,-1,1.8],[0.01,0.3,0.02],[true,false,true],['available','available','available'],[true,true,true],['0007','0007','0007'],['evaluated','evaluated','evaluated']]};",
    "class Element {",
    "  constructor(){this.children=[];this.listeners={};this.value='';this.checked=false;this.hidden=true;this.textContent='';}",
    "  append(...children){this.children.push(...children);}",
    "  replaceChildren(...children){this.children=children;}",
    "  addEventListener(type,fn){this.listeners[type]=fn;}",
    "  removeAttribute(name){delete this[name];}",
    "}",
    "function load(query,denyHistory=false){",
    "  const nodes=new Map();const el=id=>{if(!nodes.has(id))nodes.set(id,new Element());return nodes.get(id);};",
    "  el('gene-evidence-data').textContent=JSON.stringify(data);",
    "  const url=new URL('file:///D:/evaluation/report_pages/gene_evidence/index.html?'+query+'#gene-results');",
    "  const location={pathname:url.pathname,search:url.search,hash:url.hash};let writes=0;",
    "  const history={replaceState(_,__,relative){if(denyHistory)throw Error('file history denied');const next=new URL(relative,url);location.pathname=next.pathname;location.search=next.search;location.hash=next.hash;writes++;}};",
    "  vm.runInNewContext(source,{document:{getElementById:el,createElement:()=>new Element(),addEventListener(){}},",
    "    window:{},location,history,URLSearchParams,URL,Blob,setTimeout});",
    "  return {el,location,query:()=>new URLSearchParams(location.search),writes:()=>writes,",
    "    change(id,value){if(typeof value==='boolean')el(id).checked=value;else el(id).value=value;el(id).listeners.change();}};",
    "}",
    "function originQuery(side){return new URLSearchParams({gene:'001+A&B',tier:'core',analysis_id:side==='a'?A:B,collection:'GOBP-C2',",
    "  return_contrast_scope:'C01_custom output',return_contrast_id:'C01',return_category_id:'0007',return_analysis_a:A,return_analysis_b:B,return_side:side});}",
    "function returnTarget(ui){assert.equal(ui.el('contrast-return').hidden,false);return new URL(ui.el('contrast-return-link').href,'file:///D:/evaluation/report_pages/gene_evidence/index.html');}",
    "for(const side of ['a','b']){",
    "  const legacy=originQuery(side),ui=load(legacy);const initial=returnTarget(ui).href;",
    "  assert.equal(ui.el('gene-title').textContent,'001+A&B');assert.equal(ui.el('analysis-filter').value,side==='a'?A:B);",
    "  assert.equal(ui.query().get('return_analysis_id'),side==='a'?A:B);assert.equal(ui.query().get('return_collection'),'GOBP-C2');assert.equal(ui.query().get('return_tier'),'core');",
    "  ui.change('analysis-filter',C);assert.equal(ui.el('collection-filter').value,'');",
    "  ui.change('collection-filter','GOCC');ui.change('include-all',true);ui.change('leading-only',true);",
    "  assert.equal(ui.query().get('analysis_id'),C);assert.equal(ui.query().get('collection'),'GOCC');assert.equal(ui.query().get('include_all'),'1');assert.equal(ui.query().get('leading_only'),'1');",
    "  assert.equal(returnTarget(ui).href,initial);assert.equal(ui.location.hash,'#gene-results');",
    "  const reloaded=load(ui.query());assert.equal(reloaded.el('scope-error').hidden,true);assert.equal(reloaded.el('gene-title').textContent,'001+A&B');",
    "  assert.equal(reloaded.el('analysis-filter').value,C);assert.equal(reloaded.el('collection-filter').value,'GOCC');",
    "  assert.equal(reloaded.el('include-all').checked,true);assert.equal(reloaded.el('leading-only').checked,true);",
    "  assert.equal(reloaded.el('page-state').textContent,ui.el('page-state').textContent);assert.equal(returnTarget(reloaded).href,initial);",
    "  reloaded.change('collection-filter','GOBP-C2');assert.equal(reloaded.el('analysis-filter').value,'');",
    "  reloaded.change('collection-filter','');reloaded.change('include-all',false);reloaded.change('leading-only',false);",
    "  assert.equal(reloaded.query().has('analysis_id'),false);assert.equal(reloaded.query().has('collection'),false);",
    "  assert.equal(reloaded.query().has('include_all'),false);assert.equal(reloaded.query().has('leading_only'),false);",
    "  const all=load(reloaded.query());assert.equal(all.el('analysis-filter').value,'');assert.equal(all.el('collection-filter').value,'');assert.equal(returnTarget(all).href,initial);",
    "  for(const [key,value]of[['return_analysis_id',C],['return_collection','GOCC'],['return_tier','expanded']]){",
    "    const bad=new URLSearchParams(ui.query());bad.set(key,value);const rejected=load(bad);assert.equal(rejected.el('scope-error').hidden,false);assert.equal(rejected.el('contrast-return').hidden,true);assert.equal(rejected.writes(),0);",
    "  }",
    "  const partial=new URLSearchParams(legacy);partial.set('return_collection','GOBP-C2');assert.equal(load(partial).el('scope-error').hidden,false);",
    "}",
    "const independent=load(new URLSearchParams({gene:'001+A&B',include_all:'true',leading_only:'true'}));",
    "assert.equal(independent.el('contrast-return').hidden,true);assert.equal(independent.el('include-all').checked,true);assert.equal(independent.el('leading-only').checked,true);",
    "independent.change('collection-filter','GOCC');assert.equal(independent.query().get('collection'),'GOCC');assert.equal(independent.query().has('return_collection'),false);",
    "const manual=load(new URLSearchParams());manual.change('analysis-filter',C);manual.change('include-all',true);",
    "const manualReload=load(manual.query());assert.equal(manualReload.el('analysis-filter').value,C);assert.equal(manualReload.el('include-all').checked,true);assert.equal(manualReload.query().has('gene'),false);",
    "const blocked=load(originQuery('a'),true);blocked.change('analysis-filter',C);blocked.change('collection-filter','GOCC');assert.equal(blocked.el('gene-title').textContent,'001+A&B');assert.equal(blocked.writes(),0);assert.equal(returnTarget(blocked).hash.includes('category_id=0007'),true);",
    "process.stdout.write('PASS filter events, reload identity, independent A/B origins, defaults and file-history fallback\\n');"
  ), probe, useBytes = TRUE)
  result <- suppressWarnings(system2(node, c(shQuote(probe), shQuote(script)), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS filter events, reload identity, independent A/B origins, defaults and file-history fallback", fixed = TRUE)
})
