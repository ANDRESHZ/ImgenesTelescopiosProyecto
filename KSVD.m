function fk = KSVD(g,k,sigma,U,V)
    fk = zeros(size(g));
    for i = 1:k
        fk = fk + 1/sigma(i)*(U(:,i)'.*g).*V(:,i);
    end
end