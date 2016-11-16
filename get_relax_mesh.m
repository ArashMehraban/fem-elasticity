function relaxed_conn = get_relax_mesh(msh)
     
    elem_type = msh.num_nodes_per_elem;
    conn = msh.conn;
    i=1:size(conn,1);
    if(elem_type == 9)
        qt1 = [1 2 4 5];
        qt2 = [2 3 5 6];
        qt3 = [4 5 7 8];
        qt4 = [5 6 8 9];
        
        conn1 = conn(i,qt1);
        conn2 = conn(i,qt2);
        conn3 = conn(i,qt3);
        conn4 = conn(i,qt4);
        relaxed_conn = [conn1;conn2;conn3;conn4];       
    end
    if(elem_type == 27)
        cb1 = [1 2 4 5 10 11 13 14];
        cb2 = [2 3 5 6 11 12 14 15];
        cb3 = [4 5 7 8 13 14 16 17];
        cb4 = [5 6 8 9 14 15 17 18];
        cb5 = [10 11 13 14 19 20 22 23];
        cb6 = [11 12 14 15 20 21 23 24];
        cb7 = [13 14 16 17 22 23 25 26];
        cb8 = [14 15 17 18 23 24 26 27];
        
        conn1 = conn(i,cb1);
        conn2 = conn(i,cb2);
        conn3 = conn(i,cb3);
        conn4 = conn(i,cb4);
        conn5 = conn(i,cb5);
        conn6 = conn(i,cb6);
        conn7 = conn(i,cb7);
        conn8 = conn(i,cb8);
        relaxed_conn = [conn1;conn2;conn3;conn4;conn5;conn6;conn7;conn8];        
    end

end