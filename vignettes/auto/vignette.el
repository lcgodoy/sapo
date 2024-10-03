;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "vignette"
 (lambda ()
   (LaTeX-add-bibitems
    "godoy2022testing"
    "mcgarigal2013landscape"))
 '(or :bibtex :latex))

