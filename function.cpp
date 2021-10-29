#include "DVAP_par.hpp"

void adjacency_list( Instance_data &instance, Graph &G ){
	

	try{
	
	

		int    i, j;
		int    nb(0);
		
		G.lp.at( 0 ) = 0;
		
		if (Graph::oriented == 1) {
			for (i = 0; i < Graph::n; i++) {
				for (j = 0; j < Graph::m; j++) {
					if (Graph::orr.at( j ) == i) {
					G.la.at( nb ) = j;
					nb++;
					}
				}
				G.lp.at( i + 1 ) = nb;
			}
		}
		else {
			for (i = 0; i < Graph::n; i++) {
				for (j = 0; j < Graph::m; j++) {
					if (Graph::orr.at( j ) == i || Graph::ex.at( j ) == i) {
					G.la.at( nb ) = j;
					nb++;
					}
				}
				G.lp.at( i + 1 ) = nb;
			}
		}
			
			
	}
	catch (const exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		cerr << "\n\nadjacency_list - " << exc.what() << flush;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "\n\nadjacency_list - Unknown exception caught" << endl << flush;
		exit(EXIT_FAILURE);
	}


}
	

void succesor_list( Instance_data &instance, Graph &G ){

	
	try{
		

		int    i, k, u;
		
		for (k = 0; k < Graph::n; k++) {
			for (i = G.lp.at( k ); i < G.lp.at( k + 1 ); i++) {
			u = G.la.at( i );
			if ((G.ls.at( i ) = Graph::orr.at( u )) == k)
				G.ls.at( i ) = Graph::ex.at( u );
			}
		}
		
		
	}
	catch (const exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		cerr << "\n\nsuccesor_list - " << exc.what() << flush;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "\n\nsuccesor_list - Unknown exception caught" << endl << flush;
		exit(EXIT_FAILURE);
	}
}

void set_network( Instance_data &instance, Graph &G, int oriented){

	try{
		
			cout << "\nInto set_network()" << flush;
			
			int i(0),j(0),t(0), v(0), I(instance.I), J(instance.I), T(instance.T), V(instance.V);
			/* orr and ex define the arc orr-ex */
			/* orr and ex define the arc orr-ex */
			/* orr and ex define the arc orr-ex */
			

			Graph::n = (I)*(T+1);      /* nodes of the graph; a node is a pair (i,t) */
			Graph::m = (I)*(J)*(T);  /* arcs of the graph */
			Graph::I = I; Graph::T = T;Graph::V = V;
			

			Graph::oriented = oriented;
			Graph::orr = vec1Int( Graph::m);		
			Graph::ex =	 vec1Int( Graph::m); 
			

		
			/* the next code puts in order (count) head(ex) and tail(orr) of each arc orr-ex */
	
			
			
			//G.cost=	  vector2Dou( V , vector1Dou() );
			
			
			/* the next code puts in order (count) head(ex) and tail(orr) of each arc orr-ex */
			int count = 0;
			for(t = 0; t < T; ++t){
				for(i = 0; i < I; ++i){
					for(j = 0; j < J; ++j){
						if ( i == j ){
							Graph::orr.at(count) = i + I * t;
							Graph::ex.at(count)  = j + (I * (t+1));
							count++;				
						}else if ( (t + instance.tau.at(i).at(j)) <= T-1 ){  
							Graph::orr.at(count) = i + I * t;
							Graph::ex.at(count)  = j + I * (t + instance.tau.at(i).at(j));
							count++;
						}
						else if ( (t + instance.tau.at(i).at(j)) > T-1 ){  
							Graph::orr.at(count) = i + I * t;
							Graph::ex.at(count)  = (Graph::n)-(I)+j;
							count++;
						}
					}
				}
			}
	
			
			
			/* these vectors has to be updated for each vehicle in each iteration */
			G.cost=	  vec1Dou( Graph::m , 0.0) ;	 
			G.la =    vec1Int( Graph::m ) ;
			G.ls =    vec1Int( Graph::m ) ;
			G.lp =    vec1Int( Graph::n + 1) ;
			
					
			
			
			
			ofstream setred ( "red.txt", ios::out );
			setred << G;
			
			

			
			
			adjacency_list( instance, G );
			succesor_list( instance, G ); 
			
			//cout << endl;
			//copy( G.lp.begin(), G.lp.end(), output_arc );cout << endl;
			//copy( G.ls.begin(), G.ls.end(), output_arc );cout << endl;
			//copy( G.la.begin(), G.la.end(), output_arc );cout << endl;
			
			
			prec_tableau = vec1Int( Graph::n, INT_MIN );
			d_tableau =    vec1Dou( Graph::n, DBL_MAX );
			//visited =    vector1Bool( Graph::n, false );
			
			cout << "\nOut from set_network()" << flush;
			
			
		
			
	}
	catch (const exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		cerr << "\n\nset_network - " << exc.what() << flush;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "\n\nset_network - Unknown exception caught" << endl << flush;
		exit(EXIT_FAILURE);
	}
    

}


void DAG_solveshortest( Instance_data &instance, Graph &G, Veh v, Vertex orr ){
	

	
	 try{
		 

		  int i, j, k, inc_node, ii, jj, I = instance.I , J = instance.I , T = instance.T , V = instance.V;
		  Arc u;
		  double z;
		  
		  

		  //cout << "\nHello from DAG_solveshortest() for vehicle " << v << endl << flush;
		  


		  
		  /* These costs were modified according to the rule of the major*/
		 vec1Dou &cost = G.cost;

		
		/*Initialize the vector of dorrtance and the precedence vector*/
		fill ( d_tableau.begin(), d_tableau.end(), DBL_MAX);
		fill ( prec_tableau.begin(), prec_tableau.end(), INT_MIN);
		
		

		d_tableau.at(orr) = 0;
		
			
		/* Evaluate nodes in topological order */
		for ( i = orr; i < Graph::n; i++){

			inc_node = i;
			
			

			/* Evaluate adjecent nodes for each node (inc_node) */
			for ( k = G.lp.at(inc_node); k < G.lp.at(inc_node + 1); k++ ){
				
				
				

				u = G.la.at(k);
				j = G.ls.at(k);
				
				ii = inc_node%(I);
				jj = j%(J);
				
				
				
				/*Update dorrtance of nodes not permanent based on current node i*/
				z = d_tableau.at(inc_node) + cost.at(u);
				if ((z < d_tableau.at(j)) && (instance.A.at(ii).at(jj).at(v) == 1) ){
					d_tableau.at(j) = z;
					prec_tableau.at(j) = u;
				}
					
				//}
			}
			
			
		}
		//cout << endl << endl;
			
			
		//cout << "-\tBye from DAG_solveshortest()\n" << flush;

			
	}
	catch (const exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		cerr << "\n\nDAG_solveshortest - " << exc.what() << flush;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "\n\nDAG_solveshortest - Unknown exception caught" << endl << flush;
		exit(EXIT_FAILURE);
	}

	

}


