%{
     USAGE:

        chgfreq(returns, horizon , freqency, offset);

     1) offset = 0, frequency = 1:

     Takes monthly end-of month returns, in units log(R) so they may be added.
     returns monthly observations of
     k month overlapping log returns, dated as of last date.
     e.g  if r(t) is the return to end of month t, the program takes

     r(1)         0                            = ro(1)
     r(2)         0                            = ro(2)
     ..           ..
     r(k)   -->   r(1)+r(2)..+r(k)   = ro(k)
     r(k+1)       r(2)+r(3)..+r(k+1) = ro(k+1)
     ..           ..                              ..
     r(T)         r(T-k+1)+..+r(T)   = ro(T)

     NOTE returns SUMS not averages.

     2) frequency > 1, offset .ne. 0:

     samples every frequency points, starting at the frequency + offest'th.
     e.g. for freq = 3, o = 2,

     ro(1)   -->   rf(1)
     ro(2)
     ro(3)
     -----
     ro(4)   -->   rf(2)
     ro(5)
     ro(6)
       ..

    EXAMPLE:  from monthly data:
    to create quarterly data, with Q1 return = jan+feb+march,
    use horizon = 3; freqency  = 3, offset = 0;
    to create quarterly data, with Q1 return = nov+dec+jan,
    use horizon = 3; frequency = 3, offset = 2;
    to create quarterly observations of annual returns, with Q1 = feb+..+jan,
    use horizon = 12, frequency = 3, offset = 2;
    to create quarterly data, with Q1 data = march;
    use horizon = 1, frequency = 3, offset = 0;
    
    from quarterly data:
    to create annual data
    use horizon = 4, frequency = 4, offset = 0;
    to create quarterly observations of annual averages
    use horizon = 4, frequency = 1, offset = 0
 %}
    T = size(rm,1);
    ro = rm;

    if k > 1;
        bigr = rm(1:T-k+1,:);
        i = 1;
        while i <= k-1 ;
            bigr = bigr+rm(1+i:T-k+1+i,:);
            i = i+1;
        end
        ro = cat(1,(-99*ones(k-1,1))*ones(1,size(bigr,2)),bigr);
    end;

    if f > 1;
        mask = zeros(size(ro,1),1);
        i = 1;
        while i <= size(ro,1);
            if (f*i-o) <= size(ro,1);
                mask(f*i-o) = 1;
            end;
            i = i+1;
            end;
    end        
        ro = selectif(ro,mask);
    return(ro);
end