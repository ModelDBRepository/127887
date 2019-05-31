function PP = getLHSOperator( G, Np, dummyarg );

PP = G(2:Np,:) - G(1:Np-1,:);
PP = reduceByOne( PP );



