function [W,Q] = gauss_pt_wt( quadorder, qt, sdim )

% Returns the weights and coordinates of the gauss integration points

  if ( nargin < 3 )   % set default arguments
    if ( strcmp(qt,'GAUSS') == 1 )
      dim = 1;
    else
      dim = 2;
    end
  end

  if ( nargin < 2 )
    type = 'GAUSS';
  end

  if ( strcmp(qt,'GAUSS') == 1 ) 

   quadpoint=zeros(quadorder^sdim ,sdim);
   quadweight=zeros(quadorder^sdim,1);
  
    r1pt=zeros(quadorder,1); r1wt=zeros(quadorder,1);

    switch ( quadorder ) 
      case 1
        r1pt(1) = 0.000000000000000;
        r1wt(1) = 2.000000000000000;

      case 2
        r1pt(1) = 0.577350269189626;
        r1pt(2) =-0.577350269189626;

        r1wt(1) = 1.000000000000000; 
        r1wt(2) = 1.000000000000000;         

      case 3
        r1pt(1) = 0.774596669241483;
        r1pt(2) =-0.774596669241483;
        r1pt(3) = 0.000000000000000;

        r1wt(1) = 0.555555555555556;
        r1wt(2) = 0.555555555555556; 
        r1wt(3) = 0.888888888888889;   

      case 4
        r1pt(1) = 0.861134311594053;
        r1pt(2) =-0.861134311594053;
        r1pt(3) = 0.339981043584856;
        r1pt(4) =-0.339981043584856;

        r1wt(1) = 0.347854845137454;
        r1wt(2) = 0.347854845137454; 
        r1wt(3) = 0.652145154862546;
        r1wt(4) = 0.652145154862546;  

      case 5
        r1pt(1) = 0.906179845938664;
        r1pt(2) =-0.906179845938664;
        r1pt(3) = 0.538469310105683;
        r1pt(4) =-0.538469310105683;
        r1pt(5) = 0.000000000000000;

        r1wt(1) = 0.236926885056189;
        r1wt(2) = 0.236926885056189;
        r1wt(3) = 0.478628670499366;
        r1wt(4) = 0.478628670499366;  
        r1wt(5) = 0.568888888888889;  

      case 6
        r1pt(1) = 0.932469514203152;
        r1pt(2) =-0.932469514203152;
        r1pt(3) = 0.661209386466265;
        r1pt(4) =-0.661209386466265;
        r1pt(5) = 0.238619186003152;
        r1pt(6) =-0.238619186003152;

        r1wt(1) = 0.171324492379170;
        r1wt(2) = 0.171324492379170;
        r1wt(3) = 0.360761573048139;
        r1wt(4) = 0.360761573048139;   
        r1wt(5) = 0.467913934572691; 
        r1wt(6) = 0.467913934572691;
	
      case 7
        r1pt(1) =  0.949107912342759;
        r1pt(2) = -0.949107912342759;
        r1pt(3) =  0.741531185599394;
        r1pt(4) = -0.741531185599394;
        r1pt(5) =  0.405845151377397;
        r1pt(6) = -0.405845151377397;
        r1pt(7) =  0.000000000000000;

        r1wt(1) = 0.129484966168870;
        r1wt(2) = 0.129484966168870;
        r1wt(3) = 0.279705391489277;
        r1wt(4) = 0.279705391489277;
        r1wt(5) = 0.381830050505119;
        r1wt(6) = 0.381830050505119;
        r1wt(7) = 0.417959183673469;
            
       otherwise
            disp('unsupported integration order')
    
    end  % end of quadorder switch

    n=1;
     
    if ( sdim == 1 ) 
      for i = 1:quadorder
        quadpoint(n,:) = [ r1pt(i) ];           
        quadweight(n) = r1wt(i); 
        n = n+1;
      end
    
    elseif ( sdim == 2 ) 
      for i = 1:quadorder
        for j = 1:quadorder
          quadpoint(n,:) = [ r1pt(i), r1pt(j)];           
          quadweight(n) = r1wt(i)*r1wt(j); 
          n = n+1;
        end
      end  
   
    end
    
    Q=quadpoint;
    W=quadweight;
  % END OF GAUSSIAN QUADRATURE DEFINITION FOR RECTANGULAR ELEMENTS
  
  elseif ( strcmp(qt,'TRIANGULAR') == 1 )   
   
      
      if ( quadorder > 7 ) % check for valid quadrature order
        disp('Quadrature order too high for triangular quadrature');
        quadorder = 1;
      end
      
      if ( quadorder == 1 )   % set quad points and quadweights
        quadpoint = [ 0.3333333333333, 0.3333333333333 ];
        quadweight = 1;
        
      elseif ( quadorder == 2 ) 
        quadpoint = zeros( 3, 2 );
        quadweight = zeros( 3, 1 );
        
        quadpoint(1,:) = [ 0.1666666666667, 0.1666666666667 ];
        quadpoint(2,:) = [ 0.6666666666667, 0.1666666666667 ];
        quadpoint(3,:) = [ 0.1666666666667, 0.6666666666667 ]; 
        
        quadweight(1) = 0.3333333333333; 
        quadweight(2) = 0.3333333333333; 
        quadweight(3) = 0.3333333333333;   
        
      elseif ( quadorder <= 5 ) 
        quadpoint = zeros( 7, 2 );
        quadweight = zeros( 7, 1 );
        
        quadpoint(1,:) = [ 0.1012865073235, 0.1012865073235 ];
        quadpoint(2,:) = [ 0.7974269853531, 0.1012865073235 ];
        quadpoint(3,:) = [ 0.1012865073235, 0.7974269853531 ]; 
        quadpoint(4,:) = [ 0.4701420641051, 0.0597158717898 ];
        quadpoint(5,:) = [ 0.4701420641051, 0.4701420641051 ];
        quadpoint(6,:) = [ 0.0597158717898, 0.4701420641051 ]; 
        quadpoint(7,:) = [ 0.3333333333333, 0.3333333333333 ];
        
        quadweight(1) = 0.1259391805448; 
        quadweight(2) = 0.1259391805448; 
        quadweight(3) = 0.1259391805448; 
        quadweight(4) = 0.1323941527885;
        quadweight(5) = 0.1323941527885;
        quadweight(6) = 0.1323941527885;
        quadweight(7) = 0.2250000000000;  
        
      else
        quadpoint = zeros( 13, 2 );
        quadweight = zeros( 13, 1 );
        
        quadpoint(1 ,:) = [ 0.0651301029022, 0.0651301029022 ];
        quadpoint(2 ,:) = [ 0.8697397941956, 0.0651301029022 ];
        quadpoint(3 ,:) = [ 0.0651301029022, 0.8697397941956 ];
        quadpoint(4 ,:) = [ 0.3128654960049, 0.0486903154253 ];
        quadpoint(5 ,:) = [ 0.6384441885698, 0.3128654960049 ];
        quadpoint(6 ,:) = [ 0.0486903154253, 0.6384441885698 ];
        quadpoint(7 ,:) = [ 0.6384441885698, 0.0486903154253 ];
        quadpoint(8 ,:) = [ 0.3128654960049, 0.6384441885698 ];
        quadpoint(9 ,:) = [ 0.0486903154253, 0.3128654960049 ];
        quadpoint(10,:) = [ 0.2603459660790, 0.2603459660790 ];
        quadpoint(11,:) = [ 0.4793080678419, 0.2603459660790 ];
        quadpoint(12,:) = [ 0.2603459660790, 0.4793080678419 ];
        quadpoint(13,:) = [ 0.3333333333333, 0.3333333333333 ];
        
        quadweight(1 ) = 0.0533472356088;
        quadweight(2 ) = 0.0533472356088; 
        quadweight(3 ) = 0.0533472356088;
        quadweight(4 ) = 0.0771137608903;
        quadweight(5 ) = 0.0771137608903;
        quadweight(6 ) = 0.0771137608903;
        quadweight(7 ) = 0.0771137608903;
        quadweight(8 ) = 0.0771137608903;
        quadweight(9 ) = 0.0771137608903;
        quadweight(10) = 0.1756152576332; 
        quadweight(11) = 0.1756152576332; 
        quadweight(12) = 0.1756152576332;
        quadweight(13) =-0.1495700444677; 
        
      end
      
      Q=quadpoint;
      W=quadweight/2;   
    end
      
      
        
  end  % end of function
  

