
(in-package :chemicl-test)

(5am:def-suite :chemicl)
(5am:in-suite :chemicl)

(defun near (query target &optional (epsilon 0.01))
  (< (abs (- query target)) epsilon))

;;; Simple element tests

(5am:test (element.1.name)
  (let ((element (get-element 1)))
    (5am:is (equalp (name element) "hydrogen"))))

(5am:test (element.1.id)
  (let ((element (get-element 1)))
    (5am:is (equalp (id element) "H"))))

(5am:test (element.1.atomic-number)
  (let ((element (get-element 1)))
    (5am:is (= (atomic-number element) 1))))

(5am:test (element.1.mass)
  (let ((element (get-element 1)))
    (5am:is (near (mass element) 1))))

(5am:test (element.1.group)
  (let ((element (get-element 1)))
    (5am:is (= (group element) 1))))

(5am:test (element.1.period)
  (let ((element (get-element 1)))
    (5am:is (= (period element) 1))))

(5am:test element.2
  (let ((element (get-element 2)))
    (5am:is (equalp (name element) "helium"))
    (5am:is (= (atomic-number element) 2))
    (5am:is (near (mass element) 4))))

(5am:test element.3
  (let ((element (get-element "He")))
    (5am:is (equalp (name element) "helium"))
    (5am:is (= (atomic-number element) 2))
    (5am:is (near (mass element) 4))))

(5am:test element.4
  (let ((element (get-element "U")))
    (5am:is (equalp (name element) "uranium"))
    (5am:is (= (atomic-number element) 92))
    (5am:is (near (mass element) 238 .1))
    (5am:is (near (electronegativity element) 1.38))
    (5am:is (= (period element) 7))))

;;; Simple atom tests

(5am:test atom.atomic-number
  (let ((atom (make-atom "C")))
    (5am:is (= (atomic-number atom) 6))))

(5am:test atom.id
  (let ((atom (make-atom "C")))
    (5am:is (equalp (id atom) "C"))))

(5am:test atom.mass
  (let ((atom (make-atom "C")))
    (5am:is (near (mass atom) 12 0.2))))

(5am:test atom.exact-mass
  (let ((atom (make-atom "C")))
    (5am:is (near (exact-mass atom) 12 0.0001))))

(5am:test atom.get-normal-valence
  (let ((atom (make-atom "C")))
    (5am:is (equalp (get-normal-valence atom) '(4)))))

;;; Simple molecule tests

(5am:test molecule.ethane
  (let ((mol (make-molecule :name "ethane")))
    (let ((c1 (add-atom mol "C" "C1"))
          (c2 (add-atom mol "C" "C2"))
          (h1 (add-atom mol "H" "H1"))
          (h2 (add-atom mol "H" "H2"))
          (h3 (add-atom mol "H" "H3"))
          (h4 (add-atom mol "H" "H4"))
          (h5 (add-atom mol "H" "H5"))
          (h6 (add-atom mol "H" "H6")))
      (add-bond mol c1 c2)
      (add-bond mol c1 h1)
      (add-bond mol c1 h2)
      (add-bond mol c1 h3)
      (add-bond mol c2 h4)
      (add-bond mol c2 h5)
      (add-bond mol c2 h6))
    (5am:is (= (atom-count mol) 8))
    (5am:is (near (mass mol) 30.07))
    (5am:is (equalp (name mol) "ethane"))))
