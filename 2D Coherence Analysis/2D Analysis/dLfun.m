function df = dLfun(x,N)

  df = LFunc(x,N+1) - LFunc(x,N);
