


class Simulator{
	public:
		Simulator(ParameterReader parameters):parameterS(parameters){
			MaxDim.assign(3,0.) ;
			string s;
			parameterS.GetParameter("x_max",s);
			MaxDim[0] = stod(s);
			parameterS.GetParameter("y_max",s);
			MaxDim[1] = stod(s);
			parameterS.GetParameter("z_max",s);
			MaxDim[2] = stod(s);
		    MinDim.assign(3,0.) ;
		    parameterS.GetParameter("x_min",s);
			MinDim[0] = stod(s);
			parameterS.GetParameter("y_min",s);
			MinDim[1] = stod(s);
			parameterS.GetParameter("z_min",s);
			MinDim[2] = stod(s);

			parameterS.GetParameter("epsilon",s);
			epsilon = stod(s);
			parameterS.GetParameter("sigma",s);
			sigma = stod(s);
			sigma2 = sigma*sigma;
			parameterS.GetParameter("r_cut",s);
			rCut = stod(s);
			parameterS.GetParameter("t_start",s);
		    time_start = stod(s);
		    parameterS.GetParameter("t_end",s);
		    time_end = stod(s);
		    parameterS.GetParameter("delta_t",s);
		    delta = stod(s);

		    parameterS.GetParameter("vis_space",s);
			vis_space = stod(s);

			/*forn(i,3){
				cerr<<"max"<<i<<":"<<MaxDim[i]<<"\n";
				cerr<<"min"<<i<<":"<<MinDim[i]<<"\n";
			}
			cerr<<"eps:"<<epsilon<<" sig:"<<sigma<<"\n";
			cerr<<"tB:"<<time_start<<" tE:"<<time_end<<" d:"<<delta<<"\n";
			cerr<<"rCut:"<<rCut<<"\n";*/

		}

		void read(string fileName){
			ifstream in(fileName);
			double m;
			while(in>>m){
					mass.pb(m);
					vector<double> po(3,0.);
					forn(i,3)
						in>>po[i];
					position.pb(po);

					vector<double> ve(3,0.);
					forn(i,3)
						in>>ve[i];
					velocity.pb(ve);
			}
			force.assign(velocity.size(),vector<double>(3,0.));
			forceOLd.assign(velocity.size(),vector<double>(3,0.));

			cerr<<"number of particles:"<<velocity.size()<<"\n";
		}

		void bruteSimulation(){

		}

		void cellSimulation(){
			initiate();
			//return;
			double t = time_start;
			int iteration = 1, nr = 1;
			writeFile("output/",0);
			updateForces();
			while(t<=time_end){

				
				
				updatePosition();
				//return;
				updateForces();
				updateVelocities();


				if(iteration%vis_space==0){
					writeFile("output/",nr);
					++nr;
					cerr<<"iteration:"<<iteration<<"\n";
				}/*else
				if(iteration%10==0){
					cerr<<"iteration:"<<iteration<<"\n";
				}*/				

				t+=delta;				
				++iteration;
				//if(iteration==1)
				//return;
			}
		}
	private:

		void addParticleCell(int particleId){// [x,y)
			vector<size_t> dm(3,0);
			forn(i,3){
				dm[i] = (size_t)((position[particleId][i]) / rCut);
				if(dm[i]==CellDim[i])--dm[i];

				//cerr<<dm[i]<<"|"<<position[particleId][i]<<" ";
				//if(position[particleId][i] >= CellDim[i]*rCut)
				//	dm[i] = CellDim[i] - 1;//??
				//cerr<<"dim"<<i<<":"<<dm[i]<<" po:"<<position[particleId][i]<<"\n";
				//cerr<<"b:"<<(rCut*dm[i])<<" e:"<<(rCut*(dm[i]+1))<<"\n";
			}
			int id = 0;
			forn(i,3)
				id += dm[i]*CellDimCount[i];
			domain[id].pb(particleId);
			//cerr<<"id"<<id<<" prt:"<<particleId<<"\n";
		}

		void initiate(){
			CellDim.assign(3,0) ;
			CellDimCount.assign(3,0) ;
			int cells=1;
			forn(i,3){
				CellDim[i] = max(1,(int)((MaxDim[i] - MinDim[i]) / rCut));
				CellDimCount[2-i] = cells;
				cells *= CellDim[i];
				//cerr<<"dim "<<i<<":"<<CellDimCount[2-i]<<" "<<CellDim[i]<<"\n";
			}
			domain.assign(cells, Cell());
			forn(i,velocity.size()){
				addParticleCell(i);

			}
		}

		bool leaveDomain(int id, vector<double> CellMax, vector<double> CellMin){
				bool ok = false;
				forn(i,3){
					if(position[id][i] < MinDim[i]){position[id][i] += MaxDim[i]-MinDim[i]; ok = true;}
					if(position[id][i] > MaxDim[i]){position[id][i] -= MaxDim[i]-MinDim[i]; ok = true;}

					if(position[id][i] < CellMin[i]){ok = true;}
					if(position[id][i] >= CellMax[i]){ok = true;}
				}
				return ok;
		}

		void multiplyVectorScalar(std::vector<double> &u, std::vector<double> &v, double s){
			forn(i,v.size())
				u[i] = v[i]*s;
		}

		void sumVectorScalar(std::vector<double> &u, std::vector<double> &v, double s){
			forn(i,v.size())
				u[i] = v[i]+s;
		}

		void multiplyVectorVector(std::vector<double> &u1, std::vector<double> &v, std::vector<double> u){
			forn(i,v.size())
				u1[i] = v[i]*u[i];
		}

		void sumVectorVector(std::vector<double> &u1, std::vector<double> &v,  std::vector<double> u){
			forn(i,v.size())
				u1[i]=v[i]+u[i];
		}

		double Lnorm( std::vector<double> &v){
			double rs = 0.;
			forn(i,v.size())
				rs+=v[i];
			return sqrt(rs);
		}

		vector<double> getCellMin(int id){
			vector<double> d(3,0.);
			forn(i,3)
				d[i] = (((size_t)(position[id][i] / rCut))!=CellDim[i]?((int)(position[id][i] / rCut)):
					CellDim[i]-1)*rCut;
			return d;
		}

		vector<double> getCellMax(int id){
			vector<double> d(3,0.);
			forn(i,3){
				size_t k = (size_t)(position[id][i] / rCut);
				if(k < CellDim[i])
					d[i] = (k+1)*rCut;
				else d[i] = MaxDim[i];
			}
			return d;
		}

		void updatePosition(){

			vector<double> CellMax, CellMin;
			vector<int> toBeRecalculated;

			forn(j,domain.size()){
				tobeDeleted.clear();
			//	if(domain[j].size()==0)continue;

				if(domain[j].size() > 0){
					CellMax = getCellMax(*(domain[j].begin()));
					CellMin = getCellMin(*(domain[j].begin()));
					/*forn(i,3){
						cerr<<"max"<<i<<":"<<CellMax[i]<<"\n";
						cerr<<"min"<<i<<":"<<CellMin[i]<<"\n";						
					}*/
				}
				if(domain[j].size()!=0){
					//cerr<<domain[j].size()<<" id:"<<j<<" "<<"\n";
				}
				

				for (std::list<int>::iterator it=domain[j].begin(); it != domain[j].end(); ++it){
					int i1 = *it;
					/*forn(i,3)
						cerr<<position[i1][i]<<" ";
						cerr<<" id:"<<i1<<"\n";*/
					forn(i,3)
						position[i1][i] += delta*velocity[i1][i] +
						((force[i1][i] * delta*delta) / (2.*mass[i1]));

					/*forn(i,3)
						cerr<<position[i1][i]<<" ";
						cerr<<"\n";*/
					


					if(leaveDomain(i1,CellMax,CellMin))//todo - detect when leave domain and when leave just cell
					{
						//cerr<<1<<"\n";
						//continue;
						toBeRecalculated.pb(i1);
						std::list<int>::iterator it1 = it;
						tobeDeleted.pb(it1);						
					}
				}
				forn(i,tobeDeleted.size()){
					//addParticleCell(*tobeDeleted[i]);
					domain[j].erase(tobeDeleted[i]);
				}
			}
			forn(i,toBeRecalculated.size())
				addParticleCell(toBeRecalculated[i]);
		}

		void updateForcesParticle(int idCell, int idParticle){

		}

		void updatePositionNaive(){
			forn(i1,velocity.size()){
				forn(i,3){
						position[i1][i] += delta*velocity[i1][i] +
						((force[i1][i] * delta*delta) / (2.*mass[i1]));

						if(position[i1][i] < MinDim[i]) position[i1][i] += MaxDim[i];
						if(position[i1][i] > MaxDim[i]) position[i1][i] -= MaxDim[i];
				}
			}
		}

		void updateForcesNaive(){

			vector<double> offset(3,0.);

			forn(i,velocity.size()){
				forn(k,3){
					forceOLd[i][k] = force[i][k];
					force[i][k]=0.;
				}
				forn(j,velocity.size()){
					if(i==j)continue;
					double ds = 0.;//distance(i,j,offset);
					forn(k,3)
						{
							offset[k] = position[j][k]-position[i][k];
							//cerr<<[k]<<" ";
							ds+=offset[k]*offset[k];
						}
					//cerr<<ds<<"\n";
					if(ds>rCut*rCut)continue;
					double si = (sigma*sigma)/ds;
					si=si*si*si;
					//cerr<<"si"<<si<<"\n";
					double par = (1./ds)*si*(1.-2.*si);
					forn(k,3)
						force[i][k]+=24.*epsilon* par*offset[k];
				}
				//forn(k,3)
					//force[i][k]*=24.*epsilon;
			}
		}

		void updateVelocities(){
			forn(i1,position.size()){
				forn(i,3)
					velocity[i1][i] += 
					((force[i1][i] +  forceOLd[i1][i])*delta) / (2.*mass[i1]);
			}
			/*forn(i,3){
				forn(j,3)
					cerr<<"velo:"<<velocity[i][j]<<" ";
					cerr<<"\n";
			}*/
		}

		double distance(int i1,int i2, std::vector<double> offset){
			double rs = 0;
			forn(i,3){
				//assert(offset[i]==0.);
				rs += (position[i1][i]-position[i2][i] - offset[i])*(position[i1][i]-position[i2][i]-offset[i]);//??
			}
			return sqrt(rs);
		}

		vector<double> forceCell(vector<int> &po, std::vector<double> &offset, int idParticle){
			int id = 0;
			forn(i,3){
				id += CellDimCount[i]*po[i];
				if(po[i]<0){
					cerr<<"er:"<<i<<" "<<po[i]<<"\n";
					return std::vector<double>(3,0);
				}
			}

			vector<double> rs(3,0.);
			vector<double> dss(3,0.);

			for(auto it : domain[id]){
				if(it == idParticle)continue;
				
				double ds = 0.;//distance(it, idParticle, offset);
				forn(i,3){
					//if(offset[i]!=0.)
						//cerr<<offset[i]<<" i:"<<i<<"\n";
					dss[i] = (position[it][i]-offset[i])-position[idParticle][i];
					ds += dss[i]*dss[i];
				}
				double sig = sigma2/ds;
				sig = sig*sig*sig;
				double s = (1./(ds))*sig*(1.-2.*sig);
				if(ds>rCut*rCut*rCut*rCut*rCut*rCut){
					cerr<<id<<"->";
					forn(i,3)
						cerr<<po[i]<<" ";
					cerr<<"c:"<<it<<" p:"<<idParticle<<" ds:"<<ds<<"\n";
					return std::vector<double>(3,0);
				}
				forn(i,3)
					rs[i]+= s*dss[i];//??
			}
			
			return rs;
		}

		int neighbours[4]={-1,0,1};

		void getNeighboursCell(int i1, int j1, int k1){

			size_t id = CellDimCount[0]*i1+CellDimCount[1]*j1+CellDimCount[2]*k1;
			if(id >= domain.size()){
				cerr<<"EErrooorrrr!!!!!!!"<<id<<"\n";
				cerr<<i1<<" "<<j1<<" "<<k1<<"\n";
				return;
			}
			
			vector<double> offset;
			vector<int> ids;
			for(auto it : domain[id]){

				forn(i,3)
					force[it][i] = 0;

				forn(i,3)//go through neighbours
					forn(j,3)
						forn(k,3){
							//if(j==1&&i==1&&k==1)continue;

							ids.clear();
							offset.clear();
							ids.assign(3,0);
							offset.assign(3,0.);//for cell correction -offset

							int i2 = (int)(i1 + neighbours[i]);
							int j2 = (int)(j1 + neighbours[j]);
							int k2 = (int)(k1 + neighbours[k]);

							//check if boundary
							if(i2 < 0){ i2 = CellDim[0]-1; offset[0] = MaxDim[0]-MinDim[0];}
							if(j2 < 0){ j2 = CellDim[1]-1; offset[1] = MaxDim[1]-MinDim[1];}
							if(k2 < 0){ k2 = CellDim[2]-1; offset[2] = MaxDim[2]-MinDim[2];}
							if(i2 == (int)CellDim[0]){ i2 = 0; offset[0] -= MaxDim[0]-MinDim[0];}
							if(j2 == (int)CellDim[1]){ j2 = 0; offset[1] -= MaxDim[1]-MinDim[1];}
							if(k2 == (int)CellDim[2]){ k2 = 0; offset[2] -= MaxDim[2]-MinDim[2];}
							//??							
							if(i2<0 || j2<0 || k2<0){
								cerr<<i2<<" "<<j2<<" "<<k2<<"|"<<i1<<j1<<k1<<"|"<<CellDim[1]<<"\n";
								return;
							}
							ids[0]=i2;ids[1]=j2;ids[2]=k2;
							if(ids[1]<0){
								cerr<<ids[1]<<" "<<j2<<" "<<j1<<" "<<neighbours[j]<<" "<<CellDim[1];
								cerr<<" erro\n";
								return;
							}
							sumVectorVector(force[it], force[it], forceCell(ids, offset, it));
							
						}

				multiplyVectorScalar(force[it],force[it],24.*epsilon);
				
			}
		}

		void updateForces(){
			forn(i,velocity.size())
				forn(j,3)
					forceOLd[i][j] = force[i][j];

			forn(i,CellDim[0]){
				forn(j,CellDim[1]){
					forn(k,CellDim[2]){
						getNeighboursCell(i,j,k);
					}
				}
			}

			/*forn(i,3){
				forn(j,3)
					cerr<<"forces:"<<force[i][j]<<" ";
					cerr<<"\n";}*/

		}

		

		void writeFile(string path, int nr){
			string s;
			parameterS.GetParameter("name",s);
			ofstream out(path+s+to_string(nr)+".vtk");
			out<<"# vtk DataFile Version 3.0\n";
			out<<"SiWiRVisFile\n";
			out<<"ASCII\n";
			out<<"DATASET UNSTRUCTURED_GRID\n";
     		out<<"POINTS "<<position.size()<<" DOUBLE\n";
     		forn(i,position.size())
     			{
     				forn(j,3)
     					out<<position[i][j]<<" ";
     				out<<"\n";
     			}
     		out<<"POINT_DATA "<<position.size()<<"\n";
			out<<"SCALARS mass double\n";
     		out<<"LOOKUP_TABLE default\n";
     		forn(i,position.size())
     			out<<mass[i]<<"\n";
     		out<<"VECTORS force double\n";
     		forn(i,position.size())
     			{
     				forn(j,3)
     					out<<force[i][j]<<" ";
     				out<<"\n";
     			}

     		out<<"VECTORS velocity double\n";
     		forn(i,position.size())
     			{
     				forn(j,3)
     					out<<velocity[i][j]<<" ";
     				out<<"\n";
     			}

     		out<<"\n";
			out.close();
		}

		ParameterReader parameterS;

		Vector3D force, forceOLd, velocity, position;
		Vector1D mass;

		TableCells domain;// refer to id of the above vectors
		vector<double> MaxDim, MinDim;
		vector<size_t> CellDim, CellDimCount;
		double epsilon, sigma, sigma2, rCut, time_start, time_end, delta;
		int vis_space;
		vector<list<int>::iterator> tobeDeleted;
};