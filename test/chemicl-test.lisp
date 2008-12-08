
(in-package :chemicl-test)

(defparameter *tamoxifen-smiles-string*
  "CCC(C1=CC=CC=C1)=C(C2=CC=CC=C2)C3=CC=C(OCCN(C)C)C=C3")

(defparameter *smiles-test-strings*
  `(("CC(C)=O" . "CC(C)=O")
    ("CC(=O)C" . "CC(C)=O")
    (,*tamoxifen-smiles-string*
     . ,*tamoxifen-smiles-string*)))

(map nil (lambda (x)
           (destructuring-bind (smiles . result)
               x
             (assert
              (equal
               (write-smiles-string
                (parse-smiles-string 
                 smiles))
               result))))
     *smiles-test-strings*)
