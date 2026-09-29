shell_fixture_assets <- function() list(css_href = "assets/lisa_shell.css",
  js_href = "assets/lisa_shell.js", logo_href = "assets/logo.svg",
  logo_compact_href = "assets/logo-compact.svg")

test_that("shared shell shows only explicit available routes and one current link", {
  shell <- lisaR:::.lisa_report_shell(
    routes = list(overview = "../../report_index.html", analyses = list(
      href = "../single_de.html#analysis-A-GOCC", label = "Analyses"),
      contrasts = NULL, genes = "index.html", methods = ""),
    active = "genes", assets = shell_fixture_assets())
  expect_match(shell$header, 'data-lisa-route="genes" href="index.html" aria-current="page"', fixed = TRUE)
  expect_length(regmatches(shell$header, gregexpr('aria-current="page"', shell$header, fixed = TRUE))[[1L]], 1L)
  expect_false(grepl('data-lisa-route="contrasts"', shell$header, fixed = TRUE))
  expect_false(grepl('data-lisa-route="methods"', shell$header, fixed = TRUE))
  expect_match(shell$header, 'href="../../report_index.html"', fixed = TRUE)
  expect_match(shell$header, 'aria-label="Main navigation"', fixed = TRUE)
  expect_match(shell$header, 'aria-controls="lisa-shell-navigation"', fixed = TRUE)
  expect_error(lisaR:::.lisa_report_shell(routes = list(genes = "index.html"),
    active = "analyses", assets = shell_fixture_assets()), "one supplied")
  expect_error(lisaR:::.lisa_report_shell(routes = list(unknown = "x.html"),
    assets = shell_fixture_assets()), "supported route keys")
  expect_error(lisaR:::.lisa_report_shell(routes = structure(list("a", "b"),
    names = c("analyses", "analyses")), assets = shell_fixture_assets()), "uniquely named")
})

test_that("standalone shell does not invent global report links", {
  shell <- lisaR:::.lisa_report_shell(assets = shell_fixture_assets(),
    context = list(analysis = "Treatment versus baseline", collection = "GOCC"))
  expect_match(shell$header, '<span class="lisa-shell-brand">', fixed = TRUE)
  expect_false(grepl("<nav", shell$header, fixed = TRUE))
  expect_false(grepl("report_index.html", shell$header, fixed = TRUE))
  expect_match(shell$header, 'href="#lisa-main"', fixed = TRUE)
  expect_match(shell$header, 'data-lisa-context="analysis">Treatment versus baseline', fixed = TRUE)
  expect_match(shell$header, 'data-lisa-context-field="selection" hidden', fixed = TRUE)
  expect_match(shell$header, 'class="lisa-shell-identifiers" hidden', fixed = TRUE)
  expect_false(grepl("iframe", shell$header, fixed = TRUE))
  expect_match(shell$header, 'data-lisa-explorer hidden', fixed = TRUE)
})

test_that("shared branding keeps the full logo and uses its compact source only on narrow screens", {
  shell <- lisaR:::.lisa_report_shell(assets = shell_fixture_assets())
  expect_match(shell$header, 'class="lisa-shell-logo" src="assets/logo.svg"', fixed = TRUE)
  expect_match(shell$header, '<source media="(max-width: 600px)" srcset="assets/logo-compact.svg">', fixed = TRUE)
  expect_match(shell$header, 'width="410" height="128"', fixed = TRUE)
  expect_match(shell$header, 'alt="LISA - Geneset analysis, automated annotation and biological interpretation"', fixed = TRUE)
  assets <- shell_fixture_assets(); assets$logo_compact_href <- "https://example.invalid/logo.svg"
  expect_error(lisaR:::.lisa_report_shell(assets = assets), "relative offline")
  # Older callers providing one logo remain supported without invented assets.
  assets <- shell_fixture_assets(); assets$logo_compact_href <- NULL
  expect_match(lisaR:::.lisa_report_shell(assets = assets)$header,
    'srcset="assets/logo.svg"', fixed = TRUE)
})

test_that("shared shell exposes keyboard-native contextual pickers without fabricated destinations", {
  shell <- lisaR:::.lisa_report_shell(assets = shell_fixture_assets())
  expect_match(shell$header, 'role="navigation" aria-label="Explore report results"', fixed = TRUE)
  for (field in c("context", "collection", "section")) {
    expect_match(shell$header, paste0('for="lisa-nav-', field, '"'), fixed = TRUE)
    expect_match(shell$header, paste0('id="lisa-nav-', field,
      '" data-lisa-nav-select="', field, '"></select>'), fixed = TRUE)
  }
  expect_match(shell$header, 'data-lisa-nav-location aria-live="polite"', fixed = TRUE)
  expect_false(grepl('<option', shell$header, fixed = TRUE))
})

test_that("human labels and scientific context are escaped without replacing exact IDs", {
  shell <- lisaR:::.lisa_report_shell(routes = list(methods = list(href = "../data.html?a=1&b=2",
    label = "Data & methods")), active = "methods", assets = shell_fixture_assets(),
    context = list(study = "Riaz et al.", analysis = "Responders: on-treatment versus baseline",
      selection = "A < B & C", direction = "on-treatment > baseline", cutoff = "GSEA FDR <= 0.25",
      exact_ids = list(analysis_id = "longitudinal_001", category_id = "001", gene = '<script>"x"</script>')))
  expect_match(shell$header, "Responders: on-treatment versus baseline", fixed = TRUE)
  expect_match(shell$header, "A &lt; B &amp; C", fixed = TRUE)
  expect_match(shell$header, "GSEA FDR &lt;= 0.25", fixed = TRUE)
  expect_match(shell$header, '<code>001</code>', fixed = TRUE)
  expect_match(shell$header, '<code>longitudinal_001</code>', fixed = TRUE)
  expect_match(shell$header, '&lt;script&gt;&quot;x&quot;&lt;/script&gt;', fixed = TRUE)
  expect_match(shell$header, 'href="../data.html?a=1&amp;b=2"', fixed = TRUE)
  expect_false(grepl('<details class="lisa-shell-identifiers" open', shell$header, fixed = TRUE))
  expect_error(lisaR:::.lisa_report_shell(assets = shell_fixture_assets(),
    context = list(analysis = c("A", "B"))), "one non-missing")
  expect_error(lisaR:::.lisa_report_shell(assets = shell_fixture_assets(),
    context = list(exact_ids = c("001", "002"))), "unique nonempty names")
})

test_that("portable routes retain encoded scope and reject nonlocal schemes", {
  relative <- "../../../evidence/A/GOCC/index.html?category=CAT%2FA&gene=001#matrix"
  expect_identical(lisaR:::.lisa_shell_href(relative, "route"), relative)
  for (bad in c("https://example.invalid", "javascript:alert(1)", "file:///tmp/view.html",
      "/tmp/view.html", "//example.invalid", "C:\\view.html", "..\\view.html", " view.html")) {
    expect_error(lisaR:::.lisa_report_shell(routes = list(overview = bad),
      assets = shell_fixture_assets()), "relative offline")
  }
  assets <- shell_fixture_assets(); assets$js_href <- "https://example.invalid/shell.js"
  expect_error(lisaR:::.lisa_report_shell(assets = assets), "relative offline")
})

test_that("shared assets preserve existing logo bytes and produce portable hrefs", {
  source <- system.file("report_assets", package = "lisaR")
  expect_true(nzchar(source))
  out <- tempfile("shell assets with spaces ")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  assets <- lisaR:::.lisa_copy_report_shell_assets(out, asset_dir = source)
  expect_true(all(file.exists(assets$files)))
  expect_identical(assets$logo_href, "lisa-shell/LISA_logo_A1_muted_red_S_automated_annotation_final.svg")
  expect_identical(assets$logo_compact_href, "lisa-shell/LISA_logo_C_compact_icon_muted_red_S.svg")
  expect_identical(unname(tools::md5sum(file.path(source,
    "LISA_logo_A1_muted_red_S_automated_annotation_final.svg"))),
    unname(tools::md5sum(file.path(out, assets$logo_href))))
  expect_identical(unname(tools::md5sum(file.path(source,
    "LISA_logo_C_compact_icon_muted_red_S.svg"))),
    unname(tools::md5sum(file.path(out, assets$logo_compact_href))))
  expect_identical(assets, lisaR:::.lisa_copy_report_shell_assets(out, asset_dir = source))
  expect_error(lisaR:::.lisa_copy_report_shell_assets(out, asset_dir = source,
    asset_subdir = "../shell"), "descendant directory")
  expect_error(lisaR:::.lisa_copy_report_shell_assets(out, asset_dir = tempfile()), "Missing shared shell")
})

test_that("injection preserves content and attributes and creates a single main", {
  shell <- lisaR:::.lisa_report_shell(assets = shell_fixture_assets())
  content <- '<div id="view">Value \\1 &amp; exact text</div><script>var a="\\t";</script>'
  html <- paste0('<!doctype html><html><head><title>Evidence</title></head>',
    '<body class="existing" data-scope="A">', content, '</body></html>')
  actual <- lisaR:::.lisa_inject_report_shell(html, shell)
  expect_match(actual, 'class="existing lisa-shell-page" data-scope="A"', fixed = TRUE)
  expect_match(actual, content, fixed = TRUE)
  expect_match(actual, shell$head, fixed = TRUE)
  expect_match(actual, '<main id="lisa-main" class="lisa-shell-main" tabindex="-1">', fixed = TRUE)
  expect_length(regmatches(actual, gregexpr('<main', actual, fixed = TRUE))[[1L]], 1L)
  expect_error(lisaR:::.lisa_inject_report_shell(actual, shell), "already contains")
  expect_error(lisaR:::.lisa_inject_report_shell('<html><body>x</body></html>', shell), "exactly one")
})

test_that("injection labels existing main without nesting or replacing existing IDs", {
  shell <- lisaR:::.lisa_report_shell(assets = shell_fixture_assets())
  html <- '<html><head></head><body><main class="viewer">Evidence</main></body></html>'
  actual <- lisaR:::.lisa_inject_report_shell(html, shell, add_main = FALSE)
  expect_match(actual, '<main class="viewer" id="lisa-main" tabindex="-1">', fixed = TRUE)
  expect_length(regmatches(actual, gregexpr('<main', actual, fixed = TRUE))[[1L]], 1L)
  expect_error(lisaR:::.lisa_inject_report_shell(html, shell), "already has a main")
  labelled <- '<html><head></head><body><main id="content" tabindex="0">Evidence</main></body></html>'
  expect_error(lisaR:::.lisa_inject_report_shell(labelled, shell, add_main = FALSE), "another ID")
  aligned <- lisaR:::.lisa_report_shell(assets = shell_fixture_assets(), main_id = "content")
  result <- lisaR:::.lisa_inject_report_shell(labelled, aligned, add_main = FALSE)
  expect_match(result, '<main id="content" tabindex="0">Evidence</main>', fixed = TRUE)
  expect_match(result, 'href="#content"', fixed = TRUE)
})

test_that("offline enhancement keeps scoped context and neutral responsive navigation", {
  root <- system.file("report_assets", package = "lisaR")
  js <- paste(readLines(file.path(root, "lisa_shell.js"), warn = FALSE), collapse = "\n")
  css <- paste(readLines(file.path(root, "lisa_shell.css"), warn = FALSE), collapse = "\n")
  expect_match(js, "window.lisaReportShell", fixed = TRUE)
  expect_match(js, 'document.addEventListener("lisa:context"', fixed = TRUE)
  expect_match(js, 'document.dispatchEvent(new CustomEvent("lisa:shell-ready"))', fixed = TRUE)
  expect_false(grepl("fetch(", js, fixed = TRUE))
  expect_false(grepl("innerHTML", js, fixed = TRUE))
  expect_match(js, 'term.textContent = key', fixed = TRUE)
  expect_match(js, 'link.removeAttribute("aria-current")', fixed = TRUE)
  expect_match(js, 'event.key === "Escape"', fixed = TRUE)
  expect_match(css, '[aria-current]', fixed = TRUE)
  expect_match(css, '@media (max-width: 760px)', fixed = TRUE)
  expect_match(css, '@media print', fixed = TRUE)
  expect_match(css, '.lisa-shell-skip:focus', fixed = TRUE)
  expect_match(css, '.lisa-shared-shell [hidden]', fixed = TRUE)
})

test_that("an empty result inventory leaves standalone navigation quietly unavailable", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; standalone shell behavior needs its VM")
  script <- system.file("report_assets", "lisa_shell.js", package = "lisaR")
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  # Execute the complete shipped enhancement, not an extracted branch. Minimal
  # DOM support is intentional: no choices may be created for either case.
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "const source=fs.readFileSync(process.argv[2],'utf8');",
    "function run(contexts){",
    " const errors=[],bar={hidden:true,querySelector:()=>null},attrs={};",
    " const shell={querySelector:s=>s==='[data-lisa-explorer]'?bar:null,querySelectorAll:()=>[],getAttribute:k=>attrs[k],setAttribute:(k,v)=>attrs[k]=v,getBoundingClientRect:()=>({height:85})};",
    " const script={textContent:JSON.stringify({version:1,page:'overview',contexts})};",
    " const window={getComputedStyle:()=>({position:'sticky'}),addEventListener:()=>{}};",
    " const document={readyState:'complete',querySelector:s=>s==='[data-lisa-shell]'?shell:null,getElementById:id=>id==='lisa-report-navigation'?script:null,addEventListener:()=>{},dispatchEvent:()=>{},documentElement:{style:{setProperty:()=>{}}}};",
    " vm.runInNewContext(source,{window,document,console:{error:(...x)=>errors.push(x)},CustomEvent:function(type){this.type=type;}});",
    " assert.equal(window.lisaReportShell.navigation,null);assert.equal(bar.hidden,true);return errors;",
    "}",
    "assert.equal(run([]).length,0);",
    "assert.equal(run([{id:'broken',label:'Unavailable'}]).length,1);",
    "process.stdout.write('PASS empty standalone inventory without fabricated destinations\\n');"
  ), probe, useBytes = TRUE)
  result <- suppressWarnings(system2(node, c(shQuote(probe), shQuote(script)), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS empty standalone inventory", fixed = TRUE)
})

test_that("navigating to a full-report product opens its target and only its context", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; contextual details behavior needs its VM")
  script <- system.file("report_assets", "lisa_shell.js", package = "lisaR")
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "const source=fs.readFileSync(process.argv[2],'utf8');",
    "const start=source.indexOf('    function apply(choice, navigate, historyMode, moveFocus) {');",
    "const end=source.indexOf('    Object.keys(controls)',start);assert.ok(start>=0&&end>start);",
    "const context={tagName:'DETAILS',open:false,parentElement:null},otherContext={open:true};",
    "const collection={tagName:'DETAILS',open:false,parentElement:context},otherCollection={open:true};",
    "let scrolled=0,focused=0;",
    "const target={tagName:'DETAILS',open:false,parentElement:collection,closest:s=>s==='[data-lisa-nav-context]'?context:s==='[data-lisa-nav-collection]'?collection:null,scrollIntoView:()=>scrolled++,hasAttribute:()=>false,setAttribute:()=>{},focus:()=>focused++};",
    "const sandbox={changing:false,state:null,scrolling:false,external:false,refresh:()=>{},selected:()=>({section:{id:'product'}}),targetForChoice:()=>target,document:{querySelectorAll:s=>s.includes('data-lisa-nav-context')?[context,otherContext]:[collection,otherCollection]},window:{requestAnimationFrame:fn=>fn()},updateHeight:()=>{},",
    "  startPending:(t,mf)=>{t.scrollIntoView();if(mf){if(!t.hasAttribute('tabindex'))t.setAttribute('tabindex','-1');t.focus();}},stopPending:()=>{}};",
    "const choice={context:'analysis',collection:'collection',section:'product'};sandbox.choice=choice;",
    "assert.equal(vm.runInNewContext(source.slice(start,end)+';apply(choice,true,null,true)',sandbox),true);",
    "assert.equal(target.open,true);assert.equal(context.open,true);assert.equal(collection.open,true);",
    "assert.equal(otherContext.open,false);assert.equal(otherCollection.open,false);",
    "assert.equal(scrolled,1);assert.equal(focused,1);assert.equal(sandbox.changing,false);assert.equal(sandbox.scrolling,false);",
    "process.stdout.write('PASS full product target opens in its selected context\\n');"
  ), probe, useBytes = TRUE)
  result <- suppressWarnings(system2(node, c(shQuote(probe), shQuote(script)), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS full product target opens", fixed = TRUE)
})

test_that("global gene filters never leave a false fixed analysis or collection selected", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; external filter behavior needs its VM")
  script <- system.file("report_assets", "lisa_shell.js", package = "lisaR")
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "const source=fs.readFileSync(process.argv[2],'utf8'),attrs={};",
    "function element(){return {children:[],value:'',get firstChild(){return this.children[0]},appendChild(n){this.children.push(n)},removeChild(n){this.children.splice(this.children.indexOf(n),1)},addEventListener(){}}}",
    "const controls={context:element(),collection:element(),section:element()},live={},label={};",
    "const bar={hidden:true,querySelector:s=>s==='[data-lisa-nav-location]'?live:s==='[data-lisa-nav-context-label]'?label:controls[/select=\"([^\"]+)/.exec(s)?.[1]]};",
    "const shell={querySelector:s=>s==='[data-lisa-explorer]'?bar:null,querySelectorAll:()=>[],getAttribute:k=>attrs[k],setAttribute:(k,v)=>attrs[k]=v,getBoundingClientRect:()=>({height:85})};",
    "const inventory={version:1,page:'genes',mode:'external',current:{context:'A',collection:'GOCC',section:'evidence'},contexts:[{id:'A',scientific_id:'analysis_A',label:'Analysis A',kind:'analysis',collections:[{id:'GOCC',label:'GOCC',sections:[{id:'evidence',label:'Category evidence',kind:'category-evidence',href:'../single_de.html#evidence'}]}]}]};",
    "const window={getComputedStyle:()=>({position:'sticky'}),addEventListener:()=>{}};",
    "const document={readyState:'complete',querySelector:s=>s==='[data-lisa-shell]'?shell:null,getElementById:id=>id==='lisa-report-navigation'?{textContent:JSON.stringify(inventory)}:null,addEventListener:()=>{},dispatchEvent:()=>{},createElement:element,documentElement:{style:{setProperty:()=>{}}}};",
    "vm.runInNewContext(source,{window,document,console,CustomEvent:function(type){this.type=type;}});",
    "const api=window.lisaReportShell,state=()=>JSON.parse(JSON.stringify(api.navigation.getState()));",
    "assert.deepEqual(state(),{context:'A',collection:'GOCC',section:'evidence'});",
    "api.updateContext({analysis:'Analysis A',collection:'',exact_ids:{analysis_id:'analysis_A',collection:''}});",
    "assert.deepEqual(state(),{context:'A',collection:'',section:''});assert.equal(controls.collection.value,'');assert.equal(controls.section.disabled,true);assert.equal(attrs['data-navigation-aggregate'],'true');",
    "api.updateContext({analysis:'Analysis A',collection:'GOCC',exact_ids:{analysis_id:'analysis_A',collection:'GOCC'}});",
    "assert.deepEqual(state(),{context:'A',collection:'GOCC',section:'evidence'});assert.equal(attrs['data-navigation-aggregate'],'false');",
    "api.updateContext({analysis:'All individual analyses',collection:'GOCC',exact_ids:{analysis_id:'',collection:'GOCC'}});",
    "assert.deepEqual(state(),{context:'',collection:'',section:''});assert.equal(controls.context.value,'');assert.equal(attrs['data-navigation-aggregate'],'true');",
    "api.updateContext({selection:'CD8A'});assert.deepEqual(state(),{context:'',collection:'',section:''});",
    "process.stdout.write('PASS explicit all-analysis and all-collection filters clear stale scope\\n');"
  ), probe, useBytes = TRUE)
  result <- suppressWarnings(system2(node, c(shQuote(probe), shQuote(script)), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS explicit all-analysis and all-collection filters", fixed = TRUE)
})
