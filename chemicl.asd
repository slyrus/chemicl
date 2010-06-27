
(asdf:defsystem #:chemicl
  :name "chemicl"
  :author "Cyrus Harmon <ch-lisp@bobobeach.com>"
  :version "0.0.2"
  :licence "BSD"
  :description "A library for representing chemical structures"
  :depends-on (cxml cxml-stp epigraph)
  :components
  ((:cl-source-file "package")
   (:cl-source-file "util" :depends-on (package))
   (:cl-source-file "elements" :depends-on (package))
   (:cl-source-file "chemicl" :depends-on (package util elements))
   (:cl-source-file "rings" :depends-on (package elements chemicl))
   (:cl-source-file "primes" :depends-on (package elements chemicl))
   (:cl-source-file "smiles" :depends-on (package elements chemicl primes))
   (:static-file "COPYRIGHT")
   
   (:static-file "README")
   (:static-file "make-dist" :pathname #.(make-pathname :name "make-dist" :type "sh"))
   (:module "data"
            :components ((:static-file "elementdata.xml")
                         (:static-file "isotopes.xml")))))

