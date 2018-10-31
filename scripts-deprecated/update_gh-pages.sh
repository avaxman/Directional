#!/bin/bash
# Usage: cd $LIBDIRECTIONAL; scripts/update_gh-pages.sh
set -o xtrace

HEADER="title: libdirectional
author: Amir Vaxman and others
css: tutorial/style.css
html header:   <script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>
<link rel='stylesheet' href='http://yandex.st/highlightjs/7.3/styles/default.min.css'>
<script src='http://yandex.st/highlightjs/7.3/highlight.min.js'></script>
<script>hljs.initHighlightingOnLoad();</script>

"

echo "$HEADER" \
  | cat - README.md | multimarkdown -o index.html


HEADER="title: libdirectional
author: Amir Vaxman and others
css: ../tutorial/style.css
html header:   <script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>
<link rel='stylesheet' href='http://yandex.st/highlightjs/7.3/styles/default.min.css'>
<script src='http://yandex.st/highlightjs/7.3/highlight.min.js'></script>
<script>hljs.initHighlightingOnLoad();</script>

"

multimarkdown tutorial/tutorial.md -o tutorial/tutorial.html

