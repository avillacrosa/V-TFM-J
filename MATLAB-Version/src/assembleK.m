function [Kg1, Kg2] = assembleK(X,n)
	n_nodes_elem = 8;
	dim = 3;
	n_elem = size(n,1);
    Kg1 = zeros(n_nodes_elem*dim, n_elem);
    Kg2 = zeros(n_nodes_elem*dim, n_elem);
    for e = 1:n_elem
        ne = n(e,:);
        Kg1_e = zeros(n_nodes_elem*dim, 1);
        Kg2_e = zeros(n_nodes_elem*dim, 1);
        for a = 1:n_nodes_elem
            for b = 1:n_nodes_elem
                sl_a_e = (dim*(a-1)+1):dim*a;
                sl_b_e = (dim*(b-1)+1):dim*b;

                sl_a = (dim*(ne(a)-1)+1):dim*ne(a);
                sl_b = (dim*(ne(b)-1)+1):dim*ne(b);

                Kg1_e(sl_a_e) = sl_a;
                Kg2_e(sl_b_e) = sl_b;
            end
        end
        Kg1(:, e) = Kg1_e;
        Kg2(:, e) = Kg2_e;
    end
end