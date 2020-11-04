function step = m_step_fn(t,t_switchover)

% Produces a step function.  =1 before t_switchover, =0 after. 


step = floor( 1./ (floor(t./t_switchover ) + 1) );

