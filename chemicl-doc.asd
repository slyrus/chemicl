
(asdf:operate 'asdf:load-op :ch-asdf)
(asdf:operate 'asdf:load-op :smarkup)

(defpackage #:chemicl-doc-system (:use #:cl #:asdf #:ch-asdf #:smarkup))
(in-package #:chemicl-doc-system)

#.(smarkup::enable-quote-reader-macro)

(defsystem :chemicl-doc
  :name "chemicl-doc"
  :author "Cyrus Harmon" 
  :version "0.0.1"
  :depends-on (ch-asdf smarkup)
  :components
  ((:module
    "doc"
    :components
    ((:object-from-file :chemicl-doc-sexp
                        :pathname #p"chemicl-doc.sexp")

     (:filtered-object :chemicl-doc-filtered-sexp
                       :filters (:lisp :smarkup-metadata :html-metadata)
                       :depends-on (:chemicl-doc-sexp)
                       :input-object :chemicl-doc-sexp)
   
     (:filtered-object :chemicl-doc-html-filtered-sexp
                       :filters (:html-metadata)
                       :depends-on (:chemicl-doc-filtered-sexp)
                       :input-object :chemicl-doc-filtered-sexp)

     (:object-xhtml-file :chemicl-doc-xhtml
                         :pathname #p"index.xhtml"
                         :depends-on (:chemicl-doc-filtered-sexp)
                         :input-object :chemicl-doc-filtered-sexp)

     (:object-cl-pdf-file :chemicl-doc-pdf
                          :pathname #p"chemicl.pdf"
                          :depends-on (:chemicl-doc-filtered-sexp)
                          :input-object :chemicl-doc-filtered-sexp)))))

