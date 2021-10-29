#include "function.cpp"

using namespace std;
#define Variables_disaggregate

bool Generate_Instance(int argc, char *argv[], Instance_data &instance, ofstream &resultados_OV );
double DBFTVAP_solveOriginal( Instance_data &instance,  ofstream &resultados_OV  ); 	/* instance data */

int main ( int argc, char *argv[] )
{
	
	/*	List of Arguments
	 *  - 0    : program itself
	 * 	- 1    : number of terminals
	 * 	- 2    : number of periods
	 *  - 3    : number of vehicles types
	 *  - 4    : number of vehicles
	 *  - 5    : number of demands
	 *  - 6    : number of trips prohibited per type of vehicle
	 *  - 7    : file for saving objective value
	 *  - 8    : file for saving variable values
	 *  - 9    : semilla
	 *            */
	

	Instance_data instance;
		

	ofstream resultados_OV;
	ofstream resultados_VAR;
	
	
    
    if ( (argc == 10) || (argc == 9) ){
		
		resultados_OV.open( argv[7], ios::app);
		resultados_VAR.open( argv[8], ios::out);
		
		
	}else{
		
		resultados_OV.open( "resultados_OV.txt", ios::app);
		resultados_VAR.open( "resultados_VAR.txt", ios::out);
		
	}
	
	resultados_OV << "\n" << fixed << showpoint << setprecision(3) << flush;
	resultados_VAR << "\n"<< flush;	
	
	while( !Generate_Instance(argc, argv, instance, resultados_OV) )
	{
		cerr << "\nThe input files could not be read.\n" << flush;
		return 1;
	}
	
	double z  = DBFTVAP_solveOriginal( instance, resultados_OV );
	
	if ( z > 1E-3 ){
		cout << "\nDifferent Results\n" << flush;
		ofstream record_instance( "instances_bat.dat", ios::app);
		
		record_instance << "../Instances/" << n_Glo << "dist.dat" << " "
						<< "../Instances/" << n_Glo << "cos.dat" << " "
						<< "../Instances/" << n_Glo << "pro.dat" << " "
						<< "../Instances/" << n_Glo << "dem.dat" << " "
						<< "../Instances/" << n_Glo << "m.dat" << " "
						<< "../Instances/" << n_Glo << "A.dat" << " "
						<< "../Instances/" << n_Glo << ".dat" << endl << flush;
		record_instance.close();
		
		/* Both EXIT_FAILURE and the value one indicate Unsuccessful program execution status */
		return EXIT_FAILURE;
	}else{
		/* Both EXIT_SUCCESS and the value zero indicate successful program execution status */
		return EXIT_SUCCESS;
	}
	
   
   
   
}

double DBFTVAP_solveOriginal( Instance_data &ins,  ofstream &resultados_OV  ) 	/* ins data */
{           

	cout << "\nHello from DBFTVAP_solveOriginal()\n" << endl;
	
	vec3Dou &pro     = ins.pro;
	vec3Dou &cos     = ins.cos;
	vec3Dou &cosf    = ins.cosf;
	vec2Int &tau     = ins.tau;
	vec3Int &m       = ins.m;
	vec3Int &dem     = ins.dem;
	vec3Int &A       = ins.A;
	vec3Dou &In      = ins.In;
	vec2Int &Kap     = ins.Kap;
	
	
	
	int I = ins.I , J = ins.I , T = ins.T , V = ins.V;

	char var_name[1000];  
	int i,j,k,t,v;
	IloEnv env; 												// CPLEX environment (statement)

	try {
		IloModel model(env); 										// CPLEX pointer to the  model	(creates a model)

		IloArray<IloArray<IloArray<IloNumVarArray> > >							x(env);
		IloArray<IloArray<IloArray<IloNumVarArray> > >							y(env);
		IloArray<IloArray<IloNumVarArray> > 		        					q(env);
		IloArray<IloArray<IloNumVarArray> >						            	z(env);

		IloArray<IloRangeArray>                                                  Cap(env);
		IloArray<IloArray<IloRangeArray> >                                       Flux(env);
		IloArray<IloArray<IloRangeArray> >                                       BackL(env);

		IloObjective Obj;
		IloExpr objective(env);

	
		// set variable arrays x and y
		for(i = 0; i < I ; i++){
			x.add(IloArray<IloArray<IloNumVarArray> >(env));
			y.add(IloArray<IloArray<IloNumVarArray> >(env));

			for(j = 0; j < J; j++){
				x[i].add(IloArray<IloNumVarArray>(env));
				y[i].add(IloArray<IloNumVarArray>(env));

				for (t = 0; t < T; t++){
					x[i][j].add(IloNumVarArray(env));
					y[i][j].add(IloNumVarArray(env));

					for (v = 0; v < V; v++){
						sprintf(var_name, "x(%d,%d,%d,%d)", i, j, t, v);
						x[i][j][t].add(IloNumVar(env));
						x[i][j][t][v]=IloAdd( model, IloNumVar(env, 0, IloInfinity, IloNumVar::Float, var_name) );
						sprintf(var_name, "y(%d,%d,%d,%d)", i, j, t, v);
						y[i][j][t].add(IloNumVar(env));
						y[i][j][t][v]=IloAdd( model, IloNumVar(env, 0, IloInfinity, IloNumVar::Float, var_name) );
					}
				}
			}
		}
		cout << "\nX and Y ready " << flush;
		for(i = 0; i < I ; i++){
			z.add(IloArray<IloNumVarArray>(env));
			for (t = 0; t < T; t++){
				z[i].add(IloNumVarArray(env));
				for (v = 0; v < V; v++){
					sprintf(var_name, "z(%d,%d,%d)", i, t, v);
					z[i][t].add(IloNumVar(env));
					z[i][t][v]=IloAdd( model, IloNumVar(env, 0, IloInfinity, IloNumVar::Float, var_name) );
				}
			}
		}
		cout << "\nZ ready " << flush;
		for(i = 0; i < I ; i++){
			q.add(IloArray<IloNumVarArray>(env));
			for(j = 0; j < J; j++){
				q[i].add(IloNumVarArray(env));
				for (t = 0; t < T; t++){
					sprintf(var_name, "q(%d,%d,%d)", i, j, t);
					q[i][j].add(IloNumVar(env));
					q[i][j][t]=IloAdd( model, IloNumVar(env, 0, (i==j?0:IloInfinity), IloNumVar::Float, var_name) );
				}
			}
		}
		cout << "\nQ ready " << flush;
		


		cout << "\nOK1" << flush;
		//Set the objective function (Done!)
		for (i = 0; i < I; i++)
			for (j = 0; j < J; j++)
				for (t = 0; t < T; t++)
					for (v = 0; v < V; v++)	
						if (i != j ){					
							objective += pro.at(i).at(j).at(v)*x[i][j][t][v];
							objective -= cos.at(i).at(j).at(v)*y[i][j][t][v];
						}
		cout << "\nOK2" << flush;
		for (i = 0; i < I; i++)
			for (t = 0; t < T; t++)
				for (v = 0; v < V; v++)
					objective -= cosf.at(i).at(t).at(v)*z[i][t][v];	
		cout << "\nOK3" << flush;
		for (i = 0; i < I; i++)
			for (j = 0; j < J; j++)
				for (t = 0; t < T; t++)
					objective -= In.at(i).at(j).at(t)*q[i][j][t];	
		cout << "\nOK4" << flush;
		//model.add(IloMaximize(env, objective));
        Obj = IloAdd( model, IloMaximize(env, objective));
		objective.end();
		
		cout << "\nObjective Function ready " << flush;
		
		// Constraint set 1 (flow of vehicles)
		for (i = 0; i < I; i++){
			Flux.add(IloArray<IloRangeArray>(env));

			for (t = 0; t < T; t++){
				Flux[i].add(IloRangeArray(env));

				for (v = 0; v < V; v++){
					//cout << "\ni: " << i << " t: " << t << " v: " << v << flush;
					sprintf(var_name, "Flux(%d,%d,%d)", i, t, v);
					IloExpr sum_expr(env);
					Flux[i][t].add(IloRange());
					Flux[i][t][v] = IloAdd(model,IloRange(env, m.at(i).at(t).at(v), m.at(i).at(t).at(v), var_name));
					for (j = 0; j < J; j++) sum_expr += (x[i][j][t][v] + y[i][j][t][v]);
					for (k = 0; k < J; k++) if((i != k) && ((t - tau.at(k).at(i)) >= 0) ) sum_expr -= (x[k][i][t - tau.at(k).at(i)][v] + y[k][i][t - tau.at(k).at(i)][v]);
					if (t > 0) sum_expr -= y[i][i][t-1][v];
					sum_expr -= z[i][t][v];
					Flux[i][t][v].setExpr(sum_expr);
					sum_expr.end();
				}
			}
		}
		
		cout << "\nFlow Constraints ready " << flush;
		
		for (i = 0; i < I; i++){
			Cap.add(IloRangeArray(env));

			for (t = 0; t < T; t++){
				sprintf(var_name, "Cap(%d,%d)", i, t);
				IloExpr sum_expr(env);
				Cap[i].add(IloRange());
				Cap[i][t]= IloAdd(model,IloRange(env, -IloInfinity, Kap.at(i).at(t), var_name));
				for (v = 0; v < V; v++)
					for (j = 0; j < J; j++)
						 if((i != j) && ((t - tau.at(j).at(i)) >= 0) )
							 sum_expr += x[j][i][t - tau.at(j).at(i)][v] ;
				Cap[i][t].setExpr(sum_expr);
				sum_expr.end();
			}
		}
		
		cout << "\nKapacity COnstraints ready " << flush;
		
		
		cout << "\nBacklog COnstraints ready " << flush;


		// Constraint set 2 (demand)
		for(i = 0; i < I ; i++){
			BackL.add(IloArray<IloRangeArray>(env));	

			for(j = 0; j < J; j++){
				BackL[i].add(IloRangeArray(env));

				for (t = 0; t < T; t++){
					    sprintf(var_name, "BackL(%d,%d,%d)", i, j, t);
					    IloExpr sum_expr(env);
						BackL[i][j].add(IloRange());
						BackL[i][j][t] = IloAdd(model,IloRange(env, dem.at(i).at(j).at(t), dem.at(i).at(j).at(t), var_name));
						sum_expr += q[i][j][t];
						if(t > 0) sum_expr -= q[i][j][t-1];
						for(v = 0; v < V; v++) sum_expr += x[i][j][t][v];
						BackL[i][j][t].setExpr(sum_expr);
						sum_expr.end();
					}

					
			 }
		  }


		// restriction of movement on vehicles
		for (i = 0; i < I; i++)
			for (j = 0; j < J; j++)
				for (t = 0; t < T; t++)
					for (v = 0; v < V; v++)
						if ( A.at(i).at(j).at(v) == 0 ){
							x[i][j][t][v].setBounds(0,0);
							y[i][j][t][v].setBounds(0,0);
						}
					
				
			
		for (v = 0; v < V; v++)
			for (t = 0; t < T; t++)
				for (i = 0; i < I; i++)
					x[i][i][t][v].setBounds(0,0);
				
		for (i = 0; i < I; i++)
			for (j = 0; j < J; j++)
					q[i][j][T-1].setBounds(0,0);
				
		


		

		double start_cplex = clock();

		// Optimize and parameter definition
		IloCplex cplex(env);
		cplex.extract(model);
		cplex.setParam(IloCplex::EpGap,0.0);
		
		cplex.exportModel("Model_BFTVAP.lp");
		

		cplex.solve();
		
		if( cplex.getStatus() == IloAlgorithm::Optimal){
			double Z = cplex.getObjValue();
			double stop_cplex = clock();
			double solver_runtime = (stop_cplex-start_cplex)/double(CLOCKS_PER_SEC);
			cout << "Solution status = " << cplex.getStatus() << endl;
			cout << "Z  = " << cplex.getObjValue() << endl;
			cout << "time = " << solver_runtime << endl;
			cout << endl;
			
	//      op_sol << cplex.getObjValue() << endl << flush;

			IloArray<IloArray<IloArray<IloNumArray> > >							valsx(env);
			IloArray<IloArray<IloArray<IloNumArray> > >							valsy(env);
				   IloArray<IloArray<IloNumArray> > 							valsz(env);
				   IloArray<IloArray<IloNumArray> > 							valsq(env);

			for(i = 0; i < I ; i++){
				valsx.add(IloArray<IloArray<IloNumArray> >(env));
				valsy.add(IloArray<IloArray<IloNumArray> >(env));
				for(j = 0; j < J; j++){
					valsx[i].add(IloArray<IloNumArray>(env));
					valsy[i].add(IloArray<IloNumArray>(env));
					for (t = 0; t < T; t++){
						valsx[i][j].add(IloNumArray(env));		
						valsy[i][j].add(IloNumArray(env));
						for(v = 0; v < V; v++){
							valsx[i][j][t].add(IloNum());
							valsy[i][j][t].add(IloNum());
							valsx[i][j][t][v] = cplex.getValue(x[i][j][t][v]);
							valsy[i][j][t][v] = cplex.getValue(y[i][j][t][v]);
							
						}	
					}
				}
			}
			
			for(i = 0; i < I ; i++){
				valsz.add(IloArray<IloNumArray> (env));
				for (t = 0; t < T; t++){
					valsz[i].add(IloNumArray (env));
					for(v = 0; v < V; v++){
						valsz[i][t].add(IloNum());
						valsz[i][t][v] = cplex.getValue(z[i][t][v]);
						
					}
				}
			}
			
			for(i = 0; i < I ; i++){
				valsq.add(IloArray<IloNumArray> (env));
				for(j = 0; j < J; j++){
					valsq[i].add(IloNumArray (env));
					for (t = 0; t < T; t++){
						valsq[i][j].add(IloNum());
						valsq[i][j][t] = cplex.getValue(q[i][j][t]);
						
					}
				}
			}
			
			#ifdef Variables
			
			for (v = 0; v < V; v++){
				cout << endl;
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						if( valsz[i][t][v] > EPS4 ){
							cout << "\nz[" << i << "][" << t << "][" << v << "]= " << valsz[i][t][v] << flush;
							//op_sol << "\nz[" << i << "][" << t << "][" << v << "]= " << valsz[i][t][v] << flush;
							
						}
					}
				}
			}
			
			for (t = 0; t < T; t++){
				cout << endl;
				for(i = 0; i < I ; i++){
					for(j = 0; j < J; j++){	
						if( valsq[i][j][t] > EPS4 ){
							cout << "\nq[" << i << "][" << j << "][" << t << "]= " << valsq[i][j][t] << flush;
							//op_sol << "\nz[" << i << "][" << t << "][" << v << "]= " << valsz[i][t][v] << flush;
							
						}
					}
				}
			}

			for (v = 0; v < V; v++){
				cout << endl;
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						for(j = 0; j < J; j++){	
							if( valsx[i][j][t][v] > EPS4 ){
								cout << "\nx[" << i << "][" << j << "][" << t << "][" << v << "]= " << valsx[i][j][t][v] << flush;
								//cout <<setw(15) << ((t*I) + i) << setw(10) << (( (t + tau.at(i).at(j)) *I) + j)<< endl;
								//op_sol << "\nx(" << i << "," << j << "," << t << "," << v << ")= " << setprecision(4) << valsx[i][j][t][v] << flush;
							}						
						}
					}
				}
			}
			cout << endl;
			for (v = 0; v < V; v++){
				cout << endl;
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						for(j = 0; j < J; j++){	
							if( valsy[i][j][t][v] > EPS4 ){
								cout << "\ny[" << i << "][" << j << "][" << t << "][" << v << "]= " << valsy[i][j][t][v] << flush;
								//cout << setw(15) << ((t*I) + i) << setw(10) << (( (t + tau.at(i).at(j)) *I) + j)<< endl;
								//op_sol << "\ny(" << i << "," << j << "," << t << "," << v << ")= " << setprecision(4) << valsy[i][j][t][v] << flush;
								
							}
						}
					}
				}
			}
			
			#endif
		
			#ifdef Variables_disaggregate
			{
			//double tot(0.0);
			for (v = 0; v < V; v++){
				cout << endl << endl;
				//double obj_veh( 0.0 );
				for (t = 0; t < T; t++){
					for(i = 0; i < I ; i++){
						
						if( valsz[i][t][v] > EPS4 ){
							cout << "\nz[" << i << "][" << t << "][" << v << "]= " << setw(8) << valsz[i][t][v] << setw(10) << flat_ixt(ins,t,i) << flush;
							cout << setw(20) << ins.cosf.at(i).at(t).at(v) << flush;
							cout << setw(8) << ins.cosf.at(i).at(t).at(v)*( valsz[i][t][v] - ins.m.at(i).at(t).at(v)) << flush;
							//obj_veh += ins.cosf.at(i).at(t).at(v)*( valsz[i][t][v] - ins.m.at(i).at(t).at(v));
						}
						
						
						for(j = 0; j < J; j++){	
							if( valsx[i][j][t][v] > EPS4 ){
								cout << "\nx[" << i << "][" << j << "][" << t << "][" << v << "]= "<< setw(5) << valsx[i][j][t][v] << flush;
								cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								cout << setw(10) << valsx[i][j][t][v]*ins.pro.at(i).at(j).at(v) << flush;
								cout << setw(8) << ins.pro.at(i).at(j).at(v) << flush;
							}		
							if( valsy[i][j][t][v] > EPS4 ){
								cout << "\ny[" << i << "][" << j << "][" << t << "][" << v << "]= "<< setw(5) << valsy[i][j][t][v] << flush; 
								cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								cout << setw(10) << valsy[i][j][t][v]*ins.cos.at(i).at(j).at(v) << flush;
								cout << setw(8) << ins.cos.at(i).at(j).at(v) << flush;
							}
							if( valsq[i][j][t] > EPS4 ){
								cout << "\nq[" << i << "][" << j << "][" << t << "]= " << setw(5) << valsq[i][j][t] << flush;
								cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
								cout << setw(10) << valsq[i][j][t]*ins.In.at(i).at(j).at(t) << flush;
								cout << setw(8) << ins.In.at(i).at(j).at(t) << flush;
							}							
						}
					}
				}

			}

			}
			#endif
			
			
			cout << endl;
			cout << "Solution status = " << cplex.getStatus() << endl;
			cout << "Z  = " << cplex.getObjValue() << endl;
			//exit(0);

			
			resultados_OV << "For_LP:" << setw(15) << setprecision(6) << cplex.getObjValue() << "\t" << fixed << showpoint << setprecision(6) << solver_runtime <<  noshowpoint << setprecision(0) <<  flush;
			
			
		  bool flag = false;
			for (v = 0; v < V; v++){
				for(i = 0; i < I ; i++){
					for (t = 0; t < T; t++){	
						if( (valsz[i][t][v] - floor(valsz[i][t][v])) > EPS4 ) flag = true;
						for(j = 0; j < J; j++){
							if( (valsx[i][j][t][v] - floor(valsx[i][j][t][v])) > EPS4 ) flag = true;
							if( (valsy[i][j][t][v] - floor(valsy[i][j][t][v])) > EPS4 ) flag = true;
							if( (valsq[i][j][t] - floor(valsq[i][j][t])) > EPS4 ) flag = true;
						}
					}
				}
			}

			if ( flag == true){
				cout << "\n\nInstance found";
				cout << "Solution status = " << cplex.getStatus() << endl;
				cout << "\nz = " << cplex.getObjValue() << endl;
				cout << endl;
			}
	//
	//
	//		//exit(0);
	//
	//
	//		//for (i = 0; i < I; i++){
	//		//	for (j = 0; j < J; j++){
	//		//		for (t = 0; t < T; t++){				
	//		//				if(cplex.getDual(dem.at(i).at(j).at(t)) > EPS4) cout<< "DemDual[" << i+1 << "][" << j+1 << "][" << t+1 << "][" << v+1 << "] = " << cplex.getDual(dem.at(i).at(j).at(t)) << endl;
	//		//				//std::cout << "y[" << i+1 << "][" << j+1 << "][" << t+1 << "][" << v+1 << "] = " << cplex.getValue(y[i][j][t][v]) <<  std::endl;
	//		//		}
	//		//	}
	//		//}
	//		
				IloModel model2(env);
				model2.add(model);
				for(i = 0; i < I ; i++){
					for (t = 0; t < T; t++){
						for(j = 0; j < J; j++){
							model2.add(IloConversion(env, q[i][j][t], IloNumVar::Int));							
						}
						for (v = 0; v < V; v++){
							model2.add(IloConversion(env, z[i][t][v], IloNumVar::Int));	
							for(j = 0; j < J; j++){
								model2.add(IloConversion(env, x[i][j][t][v], IloNumVar::Int));							
								model2.add(IloConversion(env, y[i][j][t][v], IloNumVar::Int));
							}
						}
					}
				}
				IloCplex cplex2(env);
				cplex2.setParam(IloCplex::TiLim,10000);        
				cplex2.setParam(IloCplex::EpGap,0.0);
				cplex2.extract(model2);
				
				//// Set benders strategy to auto-generate a decomposition.
				//cplex.setParam(IloCplex::Param::Benders::Strategy,
					  //IloCplex::BendersFull);
					  
				//// Write out the auto-generated annotation.
				//cplex.writeBendersAnnotation("benders_annotations.dat");
				
				
				
				double start_cplex2 = clock();

				
				
				if( !cplex2.solve() ) {
				   env.error() << "Failed to optimize IP" << endl;
				   return -1;
				}
				
				double stop_cplex2 = clock();
				double solver_runtime2 = (stop_cplex2-start_cplex2)/double(CLOCKS_PER_SEC);
				
			#ifdef Variables
				
				for (v = 0; v < V; v++){
					cout << endl;
					for (t = 0; t < T; t++){
						for(i = 0; i < I ; i++){
							if( valsz[i][t][v] > EPS4 ){
								cout << "\nz[" << i << "][" << t << "][" << v << "]= " << valsz[i][t][v] << flush;
								//op_sol << "\nz[" << i << "][" << t << "][" << v << "]= " << valsz[i][t][v] << flush;
								
							}
						}
					}
				}
				
				for (t = 0; t < T; t++){
					cout << endl;
					for(i = 0; i < I ; i++){
						for(j = 0; j < J; j++){	
							if( valsq[i][j][t] > EPS4 ){
								cout << "\nq[" << i << "][" << j << "][" << t << "]= " << valsq[i][j][t] << flush;
								//op_sol << "\nz[" << i << "][" << t << "][" << v << "]= " << valsz[i][t][v] << flush;
								
							}
						}
					}
				}

				for (v = 0; v < V; v++){
					cout << endl;
					for (t = 0; t < T; t++){
						for(i = 0; i < I ; i++){
							for(j = 0; j < J; j++){	
								if( valsx[i][j][t][v] > EPS4 ){
									cout << "\nx[" << i << "][" << j << "][" << t << "][" << v << "]= " << valsx[i][j][t][v] << flush;
									//cout <<setw(15) << ((t*I) + i) << setw(10) << (( (t + tau.at(i).at(j)) *I) + j)<< endl;
									//op_sol << "\nx(" << i << "," << j << "," << t << "," << v << ")= " << setprecision(4) << valsx[i][j][t][v] << flush;
								}						
							}
						}
					}
				}
				cout << endl;
				for (v = 0; v < V; v++){
					cout << endl;
					for (t = 0; t < T; t++){
						for(i = 0; i < I ; i++){
							for(j = 0; j < J; j++){	
								if( valsy[i][j][t][v] > EPS4 ){
									cout << "\ny[" << i << "][" << j << "][" << t << "][" << v << "]= " << valsy[i][j][t][v] << flush;
									//cout << setw(15) << ((t*I) + i) << setw(10) << (( (t + tau.at(i).at(j)) *I) + j)<< endl;
									//op_sol << "\ny(" << i << "," << j << "," << t << "," << v << ")= " << setprecision(4) << valsy[i][j][t][v] << flush;
									
								}
							}
						}
					}
				}
				
				#endif
			
				#ifdef Variables_disaggregate
				{
				//double tot(0.0);
				for (v = 0; v < V; v++){
					cout << endl << endl;
					//double obj_veh( 0.0 );
					for (t = 0; t < T; t++){
						for(i = 0; i < I ; i++){
							
							if( valsz[i][t][v] > EPS4 ){
								cout << "\nz[" << i << "][" << t << "][" << v << "]= " << setw(8) << valsz[i][t][v] << setw(10) << flat_ixt(ins,t,i) << flush;
								cout << setw(20) << ins.cosf.at(i).at(t).at(v) << flush;
								cout << setw(8) << ins.cosf.at(i).at(t).at(v)*( valsz[i][t][v] - ins.m.at(i).at(t).at(v)) << flush;
								//obj_veh += ins.cosf.at(i).at(t).at(v)*( valsz[i][t][v] - ins.m.at(i).at(t).at(v));
							}
							
							
							for(j = 0; j < J; j++){	
								if( valsx[i][j][t][v] > EPS4 ){
									cout << "\nx[" << i << "][" << j << "][" << t << "][" << v << "]= "<< setw(5) << valsx[i][j][t][v] << flush;
									cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									cout << setw(10) << valsx[i][j][t][v]*ins.pro.at(i).at(j).at(v) << flush;
									cout << setw(8) << ins.pro.at(i).at(j).at(v) << flush;
								}		
								if( valsy[i][j][t][v] > EPS4 ){
									cout << "\ny[" << i << "][" << j << "][" << t << "][" << v << "]= "<< setw(5) << valsy[i][j][t][v] << flush; 
									cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									cout << setw(10) << valsy[i][j][t][v]*ins.cos.at(i).at(j).at(v) << flush;
									cout << setw(8) << ins.cos.at(i).at(j).at(v) << flush;
								}
								if( valsq[i][j][t] > EPS4 ){
									cout << "\nq[" << i << "][" << j << "][" << t << "]= " << setw(5) << valsq[i][j][t] << flush;
									cout << setw(10) << flat_ixt(ins,t,i) << setw(10) << flat_ixt(ins,(i==j?t+1:t + ins.tau.at(i).at(j)),j) << flush;
									cout << setw(10) << valsq[i][j][t]*ins.In.at(i).at(j).at(t) << flush;
									cout << setw(8) << ins.In.at(i).at(j).at(t) << flush;
								}							
							}
						}
					}

				}

				}
				#endif
				
			
				resultados_OV << "\tFor_IP:" << setw(15) << setprecision(6) << cplex2.getObjValue() << "\t" << fixed << showpoint << setprecision(6) << solver_runtime2 <<  noshowpoint << setprecision(0) <<  flush;
			
				
				
				double ZIP = cplex2.getObjValue();
				cout << "\nZIP  = " << ZIP << "\t"<< cplex2.getStatus() << endl;
				cout << "Z  = " << Z << "\t" << cplex.getStatus() << endl;
				double gap;
				gap = ((Z - ZIP)/ZIP)*100;
				cout << "GAP  = " << gap << " %" << endl;
				
				resultados_OV << "\t" << setw(20) << setprecision(15) << gap << flush;
				

				
				env.end();
				cout << "\nBye from FleetVap_solveOriginal()\n" << endl;
				return gap;
			}else if( cplex.getStatus() == IloAlgorithm::Infeasible || cplex.getStatus() == IloAlgorithm::Unbounded ){
				cerr << "\nProblem is infeasible or unbounded " << flush;
				cout << "\nBye from FleetVap_solveOriginal()\n" << endl;
				return 0.0;
			}
		

	}
	catch (IloException& e) {
		cerr << "\n\nConcert exception caught: " << e << endl;
		env.end();
		exit(EXIT_FAILURE);
	}
	
	catch (const exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		cerr << "\n\n" << exc.what();
		env.end();
		exit(EXIT_FAILURE);
	}

	catch (...) {
		cerr << "\n\nUnknown exception caught" << endl;
		env.end();
		exit(EXIT_FAILURE);
	}
	
		

}


bool Generate_Instance(int argc, char *argv[], Instance_data &instance, ofstream &resultados_OV ){
	
	
	cout << "\nHello from Generate_Instance\n" << flush;
	
	vec3Dou &pro     = instance.pro;
	vec3Dou &cos     = instance.cos;
	vec3Dou &cosf    = instance.cosf;
	vec2Int &tau     = instance.tau;
	vec3Int &m       = instance.m;
	vec3Int &dem     = instance.dem;
	vec3Int &A       = instance.A;
	vec3Dou &In      = instance.In;
	vec2Int &Kap     = instance.Kap;

    
  
	int terminals(25), periods(25), vehicles(25), num_vehicles(130), num_cargas(130), p(0);
	
	
	if ( (argc >= 7) ){
		terminals = atoi(argv[1]); periods = atoi(argv[2]); vehicles = atoi(argv[3]); num_vehicles = atoi(argv[4]); num_cargas = atoi(argv[5]); p = atoi(argv[6]);
	}else if ( argc == 6 ){
		terminals = atoi(argv[1]); periods = atoi(argv[2]); vehicles = atoi(argv[3]); num_vehicles = vehicles; 
		p = (terminals*terminals)/10;
	}else if(argc == 5){
		terminals = atoi(argv[1]); periods = atoi(argv[2]); vehicles = atoi(argv[3]); num_vehicles = vehicles; 
		num_cargas = (terminals*terminals*periods)/4;
		p = (terminals*terminals)/10;
	}else{
		
	}
	
	n_Glo = to_string(terminals) + "x" + to_string(periods) + "x" + to_string(vehicles) + "x" + to_string(num_vehicles) + "x" + to_string(num_cargas) + "x" + to_string(p);
	cout << "\nn_Glo: " << n_Glo << endl;
	
	string stau  ( ("../../../Documents/BFTInst/"+n_Glo + "dist.dat" ) ); 
	string  scos ( ("../../../Documents/BFTInst/"+n_Glo + "cos.dat"  ) ); 
	string  scosf( ("../../../Documents/BFTInst/"+n_Glo + "cosf.dat" ) ); 
	string spro  ( ("../../../Documents/BFTInst/"+n_Glo + "pro.dat"  ) ); 
	string sm    ( ("../../../Documents/BFTInst/"+n_Glo + "m.dat"    ) ); 
	string sdem  ( ("../../../Documents/BFTInst/"+n_Glo + "dem.dat"  ) ); 
	string sA    ( ("../../../Documents/BFTInst/"+n_Glo + "A.dat"    ) ); 
	string sH    ( ("../../../Documents/BFTInst/"+n_Glo + "H.dat"    ) ); 
	string sK    ( ("../../../Documents/BFTInst/"+n_Glo + "Kap.dat"    ) ); 
	string ss    ( ("../../../Documents/BFTInst/"+n_Glo + ".dat"     ) );
	
	cout << "\nTerminals = "		    << terminals
		 << "\nPeriods = " 		   		<< periods
		 << "\nVehicles = " 			<< vehicles 
		 << "\nNumber of Vehicles = "   << num_vehicles
		 << "\nNumber of Demands = "    << num_cargas
		 << "\nP = " 					<< p
		 << flush;
		 
		 
    if ( (terminals*terminals) < p ){
		cerr << "\nError in logic of input data, since the number of trips prohibited " << p << " can not exceed the nuumber of pair (terminalxterminal) " << (terminals*terminals) << endl;
		exit(EXIT_FAILURE);
	}
	
	if ( (terminals*terminals*periods) < num_cargas ){
		cerr << "\nError in logic of input data, since the number of freights " << num_cargas << " can not exceed the nuumber of pair (terminalxterminalxperiods) " << (terminals*terminals*periods) << endl;
		exit(EXIT_FAILURE);
	}
	if ( num_vehicles < vehicles ){
		cerr << "\nError in logic of input data, since the number of types of vehicles " << vehicles << " can not exceed the nuumber of vehicles " << num_vehicles << endl;
		exit(EXIT_FAILURE);
	}
		 
	resultados_OV << setw(3) << terminals << " " << setw(4) << periods << " " << setw(4) << vehicles << " " << setw(4) << num_vehicles << " " << setw(4) << num_cargas << " " << setw(4) << p <<"\t\t";
		 
	int I = instance.I = instance.J = terminals, T = instance.T = periods, V = instance.V = vehicles, J = terminals;
	
	if( argc == 10 ){
		semilla = atoi(argv[9]);
	}else{
		semilla = time(0);
		//semilla = 1606116100;
		
    }
	srand(semilla);
	resultados_OV << setw(10) << semilla << "\t" << flush;
	cout << "\nSemilla: " << setw(10) << semilla << flush;
	
	
	ofstream output ("Parameter_DVAP.txt", ios_base::out);
    if( !output ){
    	cerr << "File of output could not be created" << endl;
    	return false;
    }
    
    
    cout << "\nOK1" << flush;
	/* --------------------------------------------------------------- 
	 * 						Restriction of Movements
	 * --------------------------------------------------------------- */
	
	
	
	ofstream instanceData( sA.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
    
    A = vec3Int( I, vec2Int( J, vec1Int(V , 1 ) ) );
    for ( int vv = 0; vv < vehicles; vv++  ){
		
		for ( int jj = 0; jj < p; jj++ ){
			int k;
			int i = 1 + rand()%terminals;
			instanceData << i << ";";
			do{
				k = 1 + rand()%terminals;
			}while ( i == k );
			instanceData << k << ";";
			instanceData << vv+1 << ";";
			instanceData << 0 << endl;
			
			A.at(i-1).at(k-1).at(vv) = 0;
			
		}
		
	}

    instanceData.close();
    

    cout << "\nOK2" << flush;
    /* --------------------------------------------------------------- 
	 * 					Distances
	 * --------------------------------------------------------------- */
    
    instanceData.open( stau.c_str(), ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 

	tau = vec2Int( I , vec1Int( J , 0) );
	vec1Int x_cor( I );
	vec1Int y_cor( J );
	
	
	bool flag;
    do{
		flag = true;
		for( int k = 0; k < I; k++ ){
			x_cor.at(k) = 1 + rand()%(T);
			y_cor.at(k) = 1 + rand()%(T);
			//cout << "\n" << x[k] << "   " << y[k];
		}
	
	
		for(int i = 0; i < I; i++ ){
			for(int j = i; j < J; j++ ){
				if( i == j ){
					tau.at(i).at(j) = 0;
				}else{
					double val = (x_cor.at(i)^2)+(y_cor.at(j)^2);
					tau.at(i).at(j) = sqrt(val) ;
					tau.at(j).at(i) = tau.at(i).at(j);
					//tau.at(i).at(j) = 1;
					//tau.at(j).at(i) = 1;
				}
			}
		}
		
		/* Check triangular inequality for distance matrix*/
		for (int co = 0; co < I; co++ ){
			for(int i = 0; i < I; i++ ){
				for(int j = i; j < J; j++ ){
					
					if( tau.at(i).at(j) > ( tau.at(i).at(co) + tau.at(co).at(j) ) ) flag = false;
					
				}
			}
		}
	}while( flag == false );
	
	


   for(int j = 0; j < J; ++j){
		if ( j > 0 ) {instanceData << endl; }
		for(int i = 0; i < I; ++i){	   
			if ( (i != j)  && (tau.at(i).at(j) == 0)  ){
				tau.at(i).at(j) = 1;
			}
			if( tau.at(i).at(j) < 10 ){ instanceData << tau.at(i).at(j) << "   ";}else{
				instanceData << tau.at(i).at(j) << "  ";}
		}
	}
	
	instanceData.close();
	
	cout << "\nOK3" << flush;
	/* --------------------------------------------------------------- 
	 * 					Empty vehicle Costs
	 * --------------------------------------------------------------- */
	
	instanceData.open( scos.c_str(), ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 

	
	cos = vec3Dou( I, vec2Dou( J, vec1Dou(V , 0.0 ) ) );
    double maximo(0);
    int mult( 1 + rand()%20 );
	for(int vv = 0; vv < V; vv++ ){
		for(int i = 0; i < I; i++ ){
			for(int j = i; j < J; j++ ){
				
				if ( vv == 0 ){
					cos.at(i).at(j).at(0) = tau.at(i).at(j)*mult;
					cos.at(j).at(i).at(0) = cos.at(i).at(j).at(0);
				}else{
					
					cos.at(i).at(j).at(vv) = cos.at(i).at(j).at(0);
					cos.at(j).at(i).at(vv) = cos.at(i).at(j).at(vv);
				}
				
				if ( cos.at(i).at(j).at(vv) > maximo ) maximo = cos.at(i).at(j).at(vv);
				
			}
		}
	}


	
	for(int j = 0; j < I; ++j){
		if ( j > 0 ) {instanceData << endl; }
		for(int i = 0; i < I; ++i){	  
			if( cos.at(i).at(j).at(0) < 10 ){ 
				instanceData << fixed << showpoint << setprecision(1) << cos.at(i).at(j).at(0) << "   ";
			}
			else{
				instanceData << fixed << showpoint << setprecision(1) << cos.at(i).at(j).at(0) << "  ";
			}
		}
	}
	
	
	
	
	instanceData.close();
	cout << "\nOK4" << flush;
	/* --------------------------------------------------------------- 
	 * 					Fleet Sizing Costs
	 * --------------------------------------------------------------- */
	
	instanceData.open( scosf.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	cosf = vec3Dou( I, vec2Dou( T, vec1Dou(V , 0.0 ) ) );
	for (int i = 0; i < I; i++){ 
		for (int t = 0; t < T; t++){ 
			for (int v = 0; v < V; v++) {
	
				
				instanceData << i+1 << ";";
				instanceData << t+1 << ";";
				instanceData << v+1 << ";";
				double cushto(0.0);
				for (int j = 0; j < I; j++){ cushto += instance.cos.at(i).at(j).at(v); }
				cushto = cushto/I;
				cushto = ceil(((rand() / (RAND_MAX + 1.))*2)*(cushto));
				if ( cushto < EPS2 && cushto > -EPS2 ){
					
					cushto = 1+rand()%50;//cout << "\n" << cushto<< flush;
				}
				instanceData << cushto << endl;
				
				instance.cosf.at(i).at(t).at(v) = cushto;
		
	}}}
	
	instanceData.close();
	cout << "\nOK5" << flush;
	/* --------------------------------------------------------------- 
	 * 					Backlog Costs
	 * --------------------------------------------------------------- */
	
	instanceData.open( sH.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	In = vec3Dou( I, vec2Dou( J, vec1Dou(T , 0.0 ) ) );
	for (int i = 0; i < I; i++){ 
		for (int j = 0; j < I; j++){ 
			for (int t = 0; t < T; t++){ 
	
				if( i != j ){
					instanceData << i+1 << ";";
					instanceData << j+1 << ";";
					instanceData << t+1 << ";";
					double backo( 0.5 + (rand()%3) );
					instanceData << backo << endl;
				
				    instance.In.at(i).at(j).at(t) = backo;
				}
	}}}
	
	instanceData.close();
	cout << "\nOK6" << flush;
	/* --------------------------------------------------------------- 
	 * 					Terminals Capacity
	 * --------------------------------------------------------------- */
	Kap = vec2Int( I , vec1Int( T , 0) );
	instanceData.open( sK.c_str() , ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
	
	
	for (int i = 0; i < I; i++){ 
		for (int t = 0; t < T; t++){ 

	
				
				instanceData << i+1 << ";";
				instanceData << t+1 << ";";
				int Kapo( 5 + (rand()%5) );
				
				instanceData << Kapo << endl;
				
				instance.Kap.at(i).at(t) = Kapo;
		
	}}
	
	instanceData.close();
	cout << "\nOK7" << flush;
	/* --------------------------------------------------------------- 
	 * 					Profits
	 * --------------------------------------------------------------- */
	
	instanceData.open( spro.c_str(), ios::out );
	if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 

	
   pro = vec3Dou( I, vec2Dou( I, vec1Dou(V , 0.0 ) ) );

   for(int vv = 0; vv < V; vv++ ){
		for(int i = 0; i < I; i++ ){
			for(int j = i; j < I; j++ ){
				if ( vv == 0 ){				
					if( i == j ){
						pro.at(i).at(j).at(0) = 0;
					}else{
						pro.at(i).at(j).at(0) = maximo + rand()%20;
						pro.at(j).at(i).at(0) = pro.at(i).at(j).at(0);
					}			    
				}else{
					
					pro.at(i).at(j).at(vv) = pro.at(i).at(j).at(0);
					pro.at(j).at(i).at(vv) = pro.at(i).at(j).at(vv);
				}
			}
		}
	}



	for(int j = 0; j < I; ++j){
		if ( j > 0 ) {instanceData << endl; }
		for(int i = 0; i < I; ++i){	   
			if( pro.at(i).at(j).at(0) < 10 ){ 
				instanceData << fixed << showpoint << setprecision(1) << pro.at(i).at(j).at(0) << "   ";
			}else{
				instanceData << fixed << showpoint << setprecision(1) << pro.at(i).at(j).at(0) << "  ";
			}
		}
	}
	
	
	



     
    instanceData.close();
    cout << "\nOK8" << flush;
    /* --------------------------------------------------------------- 
	 * 					Supply of Vehicles
	 * --------------------------------------------------------------- */
    
    m = vec3Int( I, vec2Int( T, vec1Int(V , 0 ) ) );
    
    instanceData.open( sm.c_str(), ios::out );
    if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
    
    
    for ( int vv = 0; vv < V; vv++  ){
		
		int i = 1 + rand()%I;
		instanceData << i << ";";
		int t = 1 + rand()%(T-1);
		instanceData << t << ";";
		instanceData << vv+1 << ";";
		instanceData << vv+1 << endl;
		
		m.at(i-1).at(t-1).at(vv)++;
		
	}    
    for ( int vv = V; vv < num_vehicles; vv++  ){
		
		int i = 1 + rand()%I;
		instanceData << i << ";";
		int t = 1 + rand()%(T-1);
		instanceData << t << ";";
		instanceData << vv+1 << ";";
		int g = 1 + rand()%V;
		instanceData << g << endl;
		
		m.at(i-1).at(t-1).at(g-1)++;
		
	}
	instanceData.close();
	cout << "\nOK9" << flush;
	/* --------------------------------------------------------------- 
	 * 					Demand for Vehicles
	 * --------------------------------------------------------------- */
    

	
	
	dem = vec3Int( I, vec2Int( J, vec1Int(T , 0 ) ) );
    
    //cout << "\OK1,1" << flush; while( getchar() != '\n' );
    instanceData.open( sdem.c_str(), ios::out );
    if( !instanceData )
    {
       cerr << "File could not be opened" << endl;
       return false;
    } 
    //cout << "\OK1,2" << flush; while( getchar() != '\n' );
    for ( int n = 0; n < num_cargas; n++  ){
		int k;
		int i = 1 + rand()%I;
		instanceData << i << ";";
		do{
			k = 1 + rand()%I;
		}while ( i == k );
		instanceData << k << ";";
		int g = 1 + rand()%(T-1);
		instanceData << g << ";";
		int h = 1 + rand()%10;
		instanceData << h << endl;
		
		//cout << "\n" << i << " " << k << " " << g << " " << h;
		
		dem.at(i-1).at(k-1).at(g-1) += h;
		
	}	
	instanceData.close();
	//cout << "\OK1,3" << flush; while( getchar() != '\n' );
	
     
   	/* --------------------------------------------------------------- 
	 * 					Data Summary
	 * --------------------------------------------------------------- */
    
    instanceData.open( ss.c_str(), ios::out );
    instanceData << terminals << "\t" << periods << "\t" << vehicles << "\n" << num_vehicles << "\t" << num_cargas << "\t" << p;
    instanceData.close();
    cout << endl << flush;
    
    output << endl << endl;
    cout << "\nOK10" << flush;
    for(int x = 0; x < I; ++x){	   
	    output << endl;
		for(int y = 0; y < I; ++y){
			output << "tau["<< x+1 << "][" << y+1 << "]=" << fixed << showpoint << setprecision(1) << tau[x][y] << "  ";
		}
    }
    
    
    cout << "\nOK11" << flush;
    for (int z = 0; z < V; z++) {
	    output << endl;
		for(int y = 0; y < I; ++y){
		output << endl;
			for(int x = 0; x < I; ++x){	   
				output << "pro["<< x+1 << "][" << y+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << pro[x][y][z] << "  ";
			}
		}
    }
    output << endl << endl;
    cout << "\nOK12" << flush;
    for (int z = 0; z < V; z++) {
	    output << endl;
		for(int y = 0; y < I; ++y){
		output << endl;
			for(int x = 0; x < I; ++x){	   
				output << "cos["<< x+1 << "][" << y+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << cos[x][y][z] << "  ";
			}
		}
    }
    output << endl << endl;
    cout << "\nOK13" << flush;
    for (int z = 0; z < V; z++) {
	    output << endl;
		for(int y = 0; y < I; ++y){
		    //output << endl;
			for (int x = 0; x < T; x++){   
				output << "\ncosf["<< y+1 << "][" << x+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << cosf[y][x][z] << "  ";
			}
		}
    }
    output << endl << endl;
    cout << "\nOK14" << flush;
    for (int z = 0; z < T; z++){   
	    output << endl;
		for(int x = 0; x < I; ++x){
		//output << endl;
			for(int y = 0; y < J; ++y){
				output << "\nh["<< x+1 << "][" << y+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << In[x][y][z] << "  ";
			}
		}
    }
    output << endl << endl;
    cout << "\nOK15" << flush;
    for(int x = 0; x < I; ++x){	
	    //output << endl;
		for (int z = 0; z < T; z++){   
			output << "\nK["<< x+1 << "][" << z+1 << "]=" << fixed << showpoint << setprecision(1) << Kap[x][z] << "  ";
		}
    }
    output << endl << endl;
    cout << "\nOK16" << flush;
    for (int z = 0; z < V; z++) {
   	 //output << endl;
		for (int y = 0; y < I; y++){ 
		//output << endl;
			for (int x = 0; x < T; x++){ 
				if ( instance.m[y][x][z] != 0 ){
				output << "\nm["<< y << "][" << x << "][" << z << "]=" 
					<< setprecision(2) << instance.m[y][x][z] << setw(10) << flat_ixt(instance,x,y) << flush;
				}
			}
		}
     }
     output << endl << endl;
     cout << "\nOK17" << flush;
     for (int z = 0; z < T; z++) {
	 	//output << endl;
	 	for (int y = 0; y < I; y++){ 
	 	//output << endl;
	 		for (int x = 0; x < J; x++){ 
	 			if ( instance.dem[y][x][z] != 0 ){
	 				output << "\ndem["<< y << "][" << x << "][" << z << "]=" 
	 					<< instance.dem[y][x][z] << setw(10) << flat_ixt(instance,z,y) << setw(10) << flat_ixt(instance,(y==x?z+1:z + instance.tau.at(y).at(x)),x) << flush;
	 			}
	 		}
	 	}
      }
      output << endl << endl;
      cout << "\nOK18" << flush;
      for (int z = 0; z < V; z++) {
	  	output << endl;
	  	for (int y = 0; y < I; y++){ 
	  		output << endl;
	  		for (int x = 0; x < J; x++){ 
	  			if ( instance.A[y][x][z] == 0 ){
	  			output << "\nA["<< y << "][" << x << "][" << z << "]=" << 
	  				setprecision(2) << instance.A[y][x][z] << endl;
	  			}
	  		}
	  	}
	  }
      output << endl << endl;

     
     cout << "\nBye from Generate_Instance\n" << flush;
     
     
     return true;
     
	
}
