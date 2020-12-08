; This lisp program generates pseudocode for a function described in
; <filename>.expr, according to the following syntax:
;
;     ( (GIVEN <formal_parameters>)
;         (= <var1> <expr1>)
;         (= <var2> <expr2>) 
;         .
;         .
;         .
;	  (= <varn> <exprn>)
;	(RETURN <return_expressions>)
;     )
;
; The pseudocode is written to <filename>.pcode.  Its format is
; defined by the routine pcode_parse in pcode.h
;
; In the procedures below, the "alist" argument is an a-list mapping each
;   formal parameter or temporary variable defined so far to its register
;   address.

(defun parse (filename)
  (progn
    (setq inp_file 
      (open
        (concatenate 'string (string-downcase (string filename)) ".expr") 
        :direction :input
      )
    )
    (setq out_file 
      (open 
        (concatenate 'string (string-downcase (string filename)) ".pcode")
      	:direction :output
      )
    )
    (parse_proc (read inp_file))
    (close inp_file)
    (close out_file)
  )
)

(defun parse_proc (proc)
  (parse_body
    (cdr proc)
    (parse_GIVEN (car proc))
  )
)

(defun parse_GIVEN (cmd)
  (if (equal (car cmd) 'GIVEN)
      (progn
        (format out_file "GIVEN ~D~%" (length (cdr cmd)))
        (parse_parms (cdr cmd) nil)
      )
      (parse_error "missing GIVEN" cmd)
  )
)

(defun parse_parms (parms alist)
  (if (null parms)
      alist
      (parse_parms
        (cdr parms)
        (add_variable (car parms) alist)
      )
  )
)

(defun add_variable (var alist)
  (if (assoc var alist)
      (parse_error "redefined variable" var)
      (cons (cons var (length alist)) alist)
  )
)

(defun parse_body (body alist)
  (if (null (cdr body))
      (parse_RETURN (car body) alist)
      (parse_body
        (cdr body)
        (parse_assign (car body) alist)
      )
  )
)

(defun parse_assign (cmd alist)
  (if (equal (car cmd) '=)
      (progn
        (parse_expr (caddr cmd) alist)
        (format out_file "STORE ~D~%" (length alist))
        (add_variable (cadr cmd) alist)
      )
      (parse_error "missing =" cmd)
  )
)
  
(defun parse_RETURN (cmd alist)
  (if (equal (car cmd) 'RETURN)
      (progn
        (parse_expr_list (cdr cmd) alist)
        (format out_file "RETURN ~D~%" (length (cdr cmd)))
        nil
      )
      (parse_error "missing RETURN" cmd)
  )
)

(defun parse_expr_list (el alist)
  (if (null el)
      nil
      (progn
        (parse_expr (car el) alist)
        (parse_expr_list (cdr el) alist)
      )
  )
)

(defun parse_expr (e alist)
  (if (atom e)
      (if (numberp e)
	  (format out_file "CONST ~D~%" e)
          (parse_atom e (assoc e alist))
      )
      (progn
        (parse_expr_list (cdr e) alist)
        (format out_file "~A~%" (translate_operator (car e) (length (cdr e))))
      )
  )
)

(defun parse_atom (a apair)
  (if (null apair)
    (parse_error "undefined atom" a)
    (format out_file "LOAD ~D~%" (cdr apair))
  )
)

(defun translate_operator (op arity)
  (cond
    ((and (eq op '+)    (equal arity 2)) '+)
    ((and (eq op '-)    (equal arity 1)) 'NEG)
    ((and (eq op '-)    (equal arity 2)) '-)
    ((and (eq op '*)    (equal arity 2)) '*)
    ((and (eq op '/)    (equal arity 1)) 'INV)
    ((and (eq op '/)    (equal arity 2)) '/)
    ((and (eq op 'SQR)  (equal arity 1)) 'SQR)
    ((and (eq op 'SQRT) (equal arity 1)) 'SQRT)
    ((and (eq op 'ABS)  (equal arity 1)) 'ABS)
    ((and (eq op 'MAX)  (equal arity 2)) 'MAX)
    ((and (eq op 'MIN)  (equal arity 2)) 'MIN)
    (t (parse_error "bad or misused operator" (cons op arity)))
  )
)

(defun parse_error (msg val)
  (progn
    (print (format nil "*** SYNTAX ERROR: ~S ~A ***~%" msg val))
    nil
  )
)

