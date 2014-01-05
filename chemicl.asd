
(asdf:defsystem #:chemicl
  :name "chemicl"
  :author "Cyrus Harmon <ch-lisp@bobobeach.com>"
  :version "0.0.2"
  :licence "BSD"
  :description "A library for representing chemical structures"
  :depends-on (cxml cxml-stp epigraph parser-combinators fset)
  :in-order-to ((asdf:test-op (asdf:load-op :chemicl-test)))
  :perform (asdf:test-op :after (op c)
                         (map nil (intern (symbol-name '#:run!) '#:5am)
                              '(:chemicl :chemicl-smiles)))
  :components
  ((:cl-source-file "package")
   (:cl-source-file "util" :depends-on (package))
   (:cl-source-file "elements" :depends-on (package))
   (:cl-source-file "chemicl" :depends-on (package util elements))
   (:cl-source-file "rings" :depends-on (package elements chemicl))
   (:cl-source-file "primes" :depends-on (package elements chemicl))
   (:cl-source-file "canonicalize" :depends-on (package elements chemicl primes))
   (:cl-source-file "smiles" :depends-on (package elements chemicl primes canonicalize))
   (:static-file "COPYRIGHT")
   (:static-file "README")
   (:static-file "make-dist" :pathname #.(make-pathname :name "make-dist" :type "sh"))
   (:module "data"
            :components ((:static-file "elementdata.xml")
                         (:static-file "isotopes.xml")))))

(cl:defpackage #:chemicl-config
  (:export #:*base-directory*))

(cl:defparameter chemicl-config::*base-directory*
  (make-pathname :name nil :type nil :defaults *load-truename*))

(asdf:defsystem :chemicl-test
  :depends-on (fiveam chemicl)
  :components
  ((:module
    :test
    :components
    ((:cl-source-file "package")
     (:cl-source-file "chemicl-test" :depends-on (package))
     (:cl-source-file "smiles-test" :depends-on (package))))))
