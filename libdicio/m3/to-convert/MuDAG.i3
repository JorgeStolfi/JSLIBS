      (****************************************************************)
      (* STATE CREATION                                               *)
      (****************************************************************)
      
      (*
        These methods may either add a new state to the 
        "DAG.T", or return an existing state with the required
        properties, if it exists.  The basic implementation will
        always create a new state; derived classes may be smarter.
                
        In any case, all existing states are preserved with their
        numbers and properties.
        
        The methods will raise "Full" iff the currently allocated
        space for the DAG is exhausted (that is, if 
        "NStates() = NAlloc()" and the state is not yet present.
        In this case the client can catch the exception
        and, for instance, invoke "Alloc" or "Crunch" before 
        caling the method again.
      *)
      
      Sink(class: Symbol): State RAISES {Full};
        (*
          Returns a sink state with the given class. 
          Amortized cost: $O(1)$ time, $O(1)$ space. *)

      Append(rest: State; r: Symbol; dest: State): State RAISES {Full};
        (*
          Returns state "s" such that "Rest(s) = rest" and "Last(s) = {r,dest}". 
          Amortizeded cost: $O(1)$ time, $O(1)$ space. *)

       Copy := NIL;  

       Set := NIL;
         (*
           Disabled; use "Append" and "Sink" instead. *)

