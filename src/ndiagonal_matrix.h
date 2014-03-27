typedef struct {
  double * p_to_begin;
  unsigned side;
  unsigned max_diagonal;
  unsigned min_diagonal;

  (void *) (*allocate) (unsigned);
  void (*deallocate) (void *);
}
