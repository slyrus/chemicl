
(in-package :chemicl-test)

(5am:def-suite :chemicl)
(5am:in-suite :chemicl)

(5am:test (chemicl.1.name)
  (let ((element (get-element 1)))
    (5am:is (equalp (name element) "hydrogen"))))

(5am:test (chemicl.1.id)
  (let ((element (get-element 1)))
    (5am:is (equalp (id element) "H"))))

(5am:test (chemicl.1.atomic-number)
  (let ((element (get-element 1)))
    (5am:is (= (atomic-number element) 1))))

(5am:test (chemicl.1.mass)
  (let ((element (get-element 1)))
    (5am:is (< (- (mass element) 1) .01))))

(5am:test (chemicl.1.group)
  (let ((element (get-element 1)))
    (5am:is (= (group element) 1))))

(5am:test (chemicl.1.period)
  (let ((element (get-element 1)))
    (5am:is (= (period element) 1))))

(5am:test chemicl.2
  (let ((element (get-element 2)))
    (5am:is (equalp (name element) "helium"))
    (5am:is (= (atomic-number element) 2))
    (5am:is (< (- (mass element) 4) .01))))

(5am:test chemicl.3
  (let ((element (get-element "He")))
    (5am:is (equalp (name element) "helium"))
    (5am:is (= (atomic-number element) 2))
    (5am:is (< (- (mass element) 4) .01))))

(5am:test chemicl.4
  (let ((element (get-element "U")))
    (5am:is (equalp (name element) "uranium"))
    (5am:is (= (atomic-number element) 92))
    (5am:is (< (- (mass element) 238) .1))
    (5am:is (< (- (electronegativity element) 1.38) .01))
    (5am:is (= (period element) 7))))
