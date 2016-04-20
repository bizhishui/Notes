package benchmark;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

import field.Field;
import field.Node;
import field.VectorField;
import inputFile.InputFile;
import mathsTools.GeometryTools;
import mesh.Elements;
import mesh.Mesh;
import mesh.Vertex;
import messenger.Messenger;
import quadrature.GaussPoints12PtsOnT1;
import quadrature.Quadrature;
import simulation.SimulationContext;

/**Soft benchmark TestRegSingleLayer.java
 * 
 * Verification the regularization source method via a single layer integration
 * with specified interface force distribution.
 * Alexander Farutin et al. 2014 JCP
 * \int G_{ij}(\mathnf{R}-\mathbf{r})f_j(\mathbf{r})d^2r=4f_i(\mathbf{R})/35
 * with f_x(\mathbf{r})=r_yr_z, f_y(\mathbf{r})=r_zr_x, f_z(\mathbf{r})=r_xr_y 
 * and G_{ij}=\frac{1}{8\pi}(\frac{\delta_{ij}}{r}+\frac{r_ir_j}{r^3}).
 * 
 * @author Jinming LYU
 * @since Apr 18, 2016
 */
public class TestRegSingleLayer {
	static SimulationContext sc = new SimulationContext();
	
	public static void main(String[] args) {
		String arg = args[0];
		if (arg.contains("-h")){
			System.out.println("Help : the only avalaible option is -r");
			System.out.println("-r : run, input file given as optionnal argument");
		}
		if (arg.contains("-r")){
			try{
				InputFile file = new InputFile(args[1]);
				try {
					file.setParametersFromFile(sc.intparams, sc.doubleparams, sc.stringparams);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}catch (ArrayIndexOutOfBoundsException e){
				System.out.println("*** No input file given, generating a new one ***");
			}

			//prepare the simulation
			Mesh m1 = prepare();

			run(m1);

			//output
			//m1.saveFile(sc.repertoire_simu.concat("champ"));
			double numErr = computeSingleLayerErr(m1);
			System.out.println("\nThe max error is "+numErr);
			writeOutMeshFields(m1);
		}

	}
	
	public static void writeOutMeshFields(Mesh m) {
		double[] pos = new double[m.getFields().get(1).getNodes().length*3];
		double[] force = new double[m.getFields().get(2).getNodes().length*3];
		double[] vel = new double[m.getFields().get(3).getNodes().length*3];
		double[] pos_lim = new double[m.getFields().get(4).getNodes().length*3];
		double[] force_lim = new double[m.getFields().get(5).getNodes().length*3];
		double[] vel_lim = new double[m.getFields().get(6).getNodes().length*3];
		for (int j=0;j<3;j++) {
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++) {
				pos[3*i+j] = m.getFields().get(1).getNodes()[i].getDofValueAt(j);
				force[3*i+j] = m.getFields().get(2).getNodes()[i].getDofValueAt(j);
				vel[3*i+j] = m.getFields().get(3).getNodes()[i].getDofValueAt(j);
				pos_lim[3*i+j] = m.getFields().get(4).getNodes()[i].getDofValueAt(j);
				force_lim[3*i+j] = m.getFields().get(5).getNodes()[i].getDofValueAt(j);
				vel_lim[3*i+j] = m.getFields().get(6).getNodes()[i].getDofValueAt(j);
			}
		}
		
		String fileOutput = m.sc.repertoire_simu.concat("meshFields_".concat(String.valueOf(m.sc.Nsub)));
		m.exportMeshAndFieldsInVtkFormat(fileOutput);
		PrintWriter write_file = null;
		try {
			write_file = new PrintWriter(m.sc.repertoire_simu.concat("geom_results_".concat(String.valueOf(m.sc.Nsub))));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		write_file
		.println("# 1: pos X \t 2: pos Y \t 3: pos Z \t 4: force X \t 5: force Y \t 6: force Z \t 7:  vel X \t "
				+ "8:  vel Y \t 9: Vel Z \t 10: pos_lim X \t 11: pos_lim Y \t 12: pos_lim Z \t 13: force_lim X"
				+ "14: force_lim Y \t 15: force_lim Z \t 16: vel_lim X \t 17: vel_lim Y \t 18: vel_lim Z");

		int length = m.getVertices().size();
		for (int i=0;i<length;i++) {
			write_file.println(pos[3*i+0]+"\t"+pos[3*i+1]+"\t"+pos[3*i+2]+"\t"+
					force[3*i+0]+"\t"+force[3*i+1]+"\t"+force[3*i+2]+"\t"+
					vel[3*i+0]+"\t"+vel[3*i+1]+"\t"+vel[3*i+2]+"\t"+
					pos_lim[3*i+0]+"\t"+pos_lim[3*i+1]+"\t"+pos_lim[3*i+2]+"\t"+
					force_lim[3*i+0]+"\t"+force_lim[3*i+1]+"\t"+force_lim[3*i+2]+"\t"+
					vel_lim[3*i+0]+"\t"+vel_lim[3*i+1]+"\t"+vel_lim[3*i+2]);
		}
		
		write_file.flush();
		write_file.close();

	}
	
	public static double computeSingleLayerErr(Mesh m){
		double err = OneErr(m,3);
		return err;
	}
	
	public static double OneErr(Mesh m, int inFieldVelocity) {
		double err = 0;
		double[] velLim = new double[m.getFields().get(0).getNodes().length*3];
		double[] velTheory = new double[m.getFields().get(0).getNodes().length*3];
		double posx,posy,posz;
		for (int i=0;i<m.getFields().get(0).getNodes().length;i++) {
			posx = m.getFields().get(4).getNodes()[i].getDofValueAt(0);
			posy = m.getFields().get(4).getNodes()[i].getDofValueAt(1);
			posz = m.getFields().get(4).getNodes()[i].getDofValueAt(2);
			velTheory[3*i+0] = 4.0/35.0*posy*posz;
			velTheory[3*i+1] = 4.0/35.0*posz*posx;
			velTheory[3*i+2] = 4.0/35.0*posx*posy;
		}
		m.computeLimitValueOfField1AndStoreInField2(inFieldVelocity, 6);
		for (int j=0;j<3;j++) {
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++) {
				velLim[3*i+j] = m.getFields().get(6).getNodes()[i].getDofValueAt(j);
			}
		}
		
		for (int i=0;i<velLim.length;i++) {
			err = Math.max(err, Math.abs(velLim[i]-velTheory[i]));
		}
		return err;
	}
	
	public static void run(Mesh m1) {
		InputFile file_resume = new InputFile(sc.repertoire_simu.concat("inputFile.txt"));
		file_resume.printParametersInFile(sc.intparams, sc.doubleparams, sc.stringparams);
		//compute the pos_lim
        //m1.computeLimitValueOfField1AndStoreInField2(1, 4);
		
		computeSurfaceVelocity(m1);
	}
	
	public static void computeSurfaceVelocity(Mesh m) {
		applyPrefixForce(m,5,2);
		computeSingleLayerVelocity(m,2,3);
		utilityCollocateVelocity(m,3,3);
	}
	
	//convert the physical space velocity in nodal values
	public static void utilityCollocateVelocity(Mesh m, final int inFieldVelocity, final int outFieldVelocity){
		
		
		double[] velocity = new double[m.getFields().get(0).getNodes().length*3];
		m.resetField(6);
		for (int j=0;j<3;j++){
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++){
				velocity[3*i+j] = m.getFields().get(inFieldVelocity).getNodes()[i].getDofValueAt(j);
				m.getFields().get(6).getNodes()[i].setDofValueAt(j, velocity[3*i+j]);
			}
		}
		
		velocity=m.convertValueAtNodeInNodalValues(velocity);
		for (int j=0;j<3;j++){
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++){
				m.getFields().get(outFieldVelocity).getNodes()[i].setDofValueAt(j, velocity[3*i+j]);
			}
		}
		
	}
	
    //compute the single layer velocity field (f3) by force field (f2), all in nodal version
	public static void computeSingleLayerVelocity(Mesh m, final int inFieldForce, final int outFieldVelocity) {
		m.resetField(outFieldVelocity);
		if (m.sc.arch.equals("CPU")) {
			if (m.sc.regularized>0) {
				System.out.println("Using Regularized Green function");
				//computeSingleLayerVelocity_CPUReg(m,inFieldForce,outFieldVelocity);
				//computeSingleLayerVelocity_1CPUReg(m,inFieldForce,outFieldVelocity);
				computeSingleLayerVelocity_1CPURegOpt(m,inFieldForce,outFieldVelocity);
			} else {
				System.out.println("Using Original Green function");
				computeSingleLayerVelocity_CPUCurrent(m,inFieldForce,outFieldVelocity);
			}
		} else {
			System.out.println("GPU for single layer computation is not ready yet!");
		}
	}
	
	//compute at each node the physical value of single-layer value
	//Only use regularized source when the node is irregular
	public static void computeSingleLayerVelocity_1CPURegOpt(Mesh m, final int inFieldForce, final int outFieldVelocity) {
		//System.out.println("This function is called!");
		final int Nnodes = m.getFields().get(inFieldForce).getNodes().length;
		final double[][] xx0 = new double[Nnodes][3];
		//compute all the nodes' coordinates (target position, R) and stored in xx0
		for (int i=0;i<Nnodes;i++) {
			Elements e = m.getFields().get(inFieldForce).getNodes()[i].getVertex().getElements().get(1);
			Node n0 = m.getFields().get(inFieldForce).getNodes()[i];
			Vertex v0 = n0.getVertex();
			double[] shape;
			
			if (e.getVertices()[0].equals(v0)){
				shape = e.computeShapeFunctionAt(new double[]{0.0,0.0});
			}else{
				if (e.getVertices()[1].equals(v0)){
					shape = e.computeShapeFunctionAt(new double[]{1.0,0.0});
				}else{
					shape = e.computeShapeFunctionAt(new double[]{0.0,1.0});
				}
			}
			//real position of vertex
			double[] x0 = new double[3];
			for (int kk=0;kk<shape.length;kk++){
				for (int jj=0;jj<x0.length;jj++){
					x0[jj]+= shape[kk]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
				}
			}
			xx0[i] = x0;
		}
		
		//Test consistency of pos_lim
//			for (int j=0;j<3;j++) {
//				for (int i=0;i<Nnodes;i++) {
//					System.out.println(Math.abs(m.getFields().get(4).getNodes()[i].getDofValueAt(j)-xx0[i][j]));
//				}
//			}
		
		Quadrature quad  = new GaussPoints12PtsOnT1();
		double G00,G01,G02,G11,G12,G22;
		double epsilon = m.sc.reg_epsilon;

		double[] shape_e;
		double[][] shape_e_der;
		double[] tksi= new double[3];
		double[] teta= new double[3];
		double[] coordinatesreal = new double[3];
		double[] f = new double[3];
		double wk;
		double[] normale;
		double g;
		double wkg;
		double gksiksi;
		double gksieta;
		double getaeta;
		Node n0;
		double[] x0;
		double[] v_n0 = new double[3];
		double[] vitn0x =new double[Nnodes];
		double[] vitn0y =new double[Nnodes];
		double[] vitn0z =new double[Nnodes];

		int nbquad = quad.getNumberOfQuadraturePoints();
		double[] y = new double[3];
		double pi8=8.0*Math.PI;
		double oneOverR;
		double oneOverR2;
		double oneOverR3;
		double R2;
		double epsCrit = Math.sqrt(2*Math.PI/5.0/Math.pow(4.0, m.sc.Nsub*1.0));

		//source position loops, r
		for (int j=0;j<m.getElements().length;j++){ 
			for (int k=0;k<nbquad;k++){

				wk = quad.getWeightPoint(k);
				shape_e = m.getElements(j).computeShapeFunctionAt(quad.getQuadraturePoint(k));
				shape_e_der = m.getElements(j).computeShapeFunctionFirstDerivativesAt(quad.getQuadraturePoint(k));
				//System.out.println("The shape func's length is "+shape_e.length);
				
				cleanDoubleArray(tksi);
				cleanDoubleArray(teta);
				cleanDoubleArray(coordinatesreal);
				cleanDoubleArray(f);

				for (int l=0;l<shape_e.length;l++){
					for (int n=0;n<coordinatesreal.length;n++){
						//m.getElements(j).getOneRing(l): return the l-th Vertex of the j-th element at mesh m
						//getNodes(1): return the 1-th node define on the vertex???
						//getDofValues(n): Returns the n-th component of the dof (nodal values) values at the node
						//coordinatesreal are the real coordinate at a Gauss quadrature point of a element, i.e., xi=Xij*Nj
						coordinatesreal[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
						//f is the real physical space force, i.e., f_lim
						f[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(inFieldForce).getDofValues(n);//inFieldForce=2

						tksi[n]+=shape_e_der[l][0]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
						teta[n]+=shape_e_der[l][1]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);

					} //for n, coordinatesreal.length
				} // for l, shape_e.length



				normale = GeometryTools.CrossProduct(tksi, teta);
				GeometryTools.Normalize(normale);

				gksiksi=GeometryTools.ScalarProduct(tksi, tksi);
				gksieta=GeometryTools.ScalarProduct(tksi, teta);
				getaeta=GeometryTools.ScalarProduct(teta, teta);

				g = gksiksi*getaeta-gksieta*gksieta;
				//wkg = wk*Math.sqrt(g);
				wkg = wk*Math.sqrt(g)/pi8;

				//target position (vertex) loops, R
				//Any boundary source will influence all the vertex in this mesh
				for (int i=0;i<Nnodes;i++){
					x0 = xx0[i];

					y[0]=coordinatesreal[0]-x0[0];
					y[1]=coordinatesreal[1]-x0[1];
					y[2]=coordinatesreal[2]-x0[2];

					R2 = y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
					int isRegular;
					if (Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])<1.5*epsCrit) {
						oneOverR =1./Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+epsilon*epsilon);
						isRegular = 0;
						//System.out.println("A node irregular.");
					} else {
						oneOverR =1./Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
						isRegular = 1;
						//System.out.println("A node regular.");
					}

					oneOverR2= oneOverR*oneOverR;
					oneOverR3= oneOverR2*oneOverR;

					if (isRegular == 1) {
						G00 = oneOverR+y[0]*y[0]*oneOverR3;
						G01 = oneOverR3*y[0]*y[1];
						G02 = oneOverR3*y[0]*y[2];
						G11 = oneOverR+y[1]*y[1]*oneOverR3;
						G12 = oneOverR3*y[1]*y[2];
						G22 = oneOverR+y[2]*y[2]*oneOverR3;
					} else {
						G00 = (R2+2*epsilon*epsilon+y[0]*y[0])*oneOverR3;
						G01 = oneOverR3*y[0]*y[1];
						G02 = oneOverR3*y[0]*y[2];
						G11 = (R2+2*epsilon*epsilon+y[1]*y[1])*oneOverR3;
						G12 = oneOverR3*y[1]*y[2];
						G22 = (R2+2*epsilon*epsilon+y[2]*y[2])*oneOverR3;
					}

					//f[] is the physical coords values of forces
					v_n0[0]= wkg*(G00*f[0]+G01*f[1]+G02*f[2]);
					v_n0[1]= wkg*(G01*f[0]+G11*f[1]+G12*f[2]);
					v_n0[2]= wkg*(G02*f[0]+G12*f[1]+G22*f[2]);


					vitn0x[i] += v_n0[0];
					vitn0y[i] += v_n0[1];
					vitn0z[i] += v_n0[2];

				}//end nodes
			}//end quad points
		}//end elements

		for (int i=0;i<Nnodes;i++){
			n0 = m.getFields(outFieldVelocity).getNodes(i);
			n0.setDofValueAt(0, n0.getDofValueAt(0)+vitn0x[i]);
			n0.setDofValueAt(1, n0.getDofValueAt(1)+vitn0y[i]);
			n0.setDofValueAt(2, n0.getDofValueAt(2)+vitn0z[i]);
		}
				
	}
		
	//compute at each node the physical value of single-layer value
	public static void computeSingleLayerVelocity_1CPUReg(Mesh m, final int inFieldForce, final int outFieldVelocity) {
		final int Nnodes = m.getFields().get(inFieldForce).getNodes().length;
		final double[][] xx0 = new double[Nnodes][3];
		//compute all the nodes' coordinates (target position, R) and stored in xx0
		for (int i=0;i<Nnodes;i++) {
			Elements e = m.getFields().get(inFieldForce).getNodes()[i].getVertex().getElements().get(1);
			Node n0 = m.getFields().get(inFieldForce).getNodes()[i];
			Vertex v0 = n0.getVertex();
			double[] shape;
			
			if (e.getVertices()[0].equals(v0)){
				shape = e.computeShapeFunctionAt(new double[]{0.0,0.0});
			}else{
				if (e.getVertices()[1].equals(v0)){
					shape = e.computeShapeFunctionAt(new double[]{1.0,0.0});
				}else{
					shape = e.computeShapeFunctionAt(new double[]{0.0,1.0});
				}
			}
			//real position of vertex
			double[] x0 = new double[3];
			for (int kk=0;kk<shape.length;kk++){
				for (int jj=0;jj<x0.length;jj++){
					x0[jj]+= shape[kk]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
				}
			}
			xx0[i] = x0;
		}
		
		//Test consistency of pos_lim
//		for (int j=0;j<3;j++) {
//			for (int i=0;i<Nnodes;i++) {
//				System.out.println(Math.abs(m.getFields().get(4).getNodes()[i].getDofValueAt(j)-xx0[i][j]));
//			}
//		}
		
		Quadrature quad  = new GaussPoints12PtsOnT1();
		double G00,G01,G02,G11,G12,G22;
		double epsilon = m.sc.reg_epsilon;

		double[] shape_e;
		double[][] shape_e_der;
		double[] tksi= new double[3];
		double[] teta= new double[3];
		double[] coordinatesreal = new double[3];
		double[] f = new double[3];
		double wk;
		double[] normale;
		double g;
		double wkg;
		double gksiksi;
		double gksieta;
		double getaeta;
		Node n0;
		double[] x0;
		double[] v_n0 = new double[3];
		double[] vitn0x =new double[Nnodes];
		double[] vitn0y =new double[Nnodes];
		double[] vitn0z =new double[Nnodes];

		int nbquad = quad.getNumberOfQuadraturePoints();
		double[] y = new double[3];
		double pi8=8.0*Math.PI;
		double oneOverR;
		double oneOverR2;
		double oneOverR3;
		double R2;

		//source position loops, r
		for (int j=0;j<m.getElements().length;j++){ 
			for (int k=0;k<nbquad;k++){

				wk = quad.getWeightPoint(k);
				shape_e = m.getElements(j).computeShapeFunctionAt(quad.getQuadraturePoint(k));
				shape_e_der = m.getElements(j).computeShapeFunctionFirstDerivativesAt(quad.getQuadraturePoint(k));

				cleanDoubleArray(tksi);
				cleanDoubleArray(teta);
				cleanDoubleArray(coordinatesreal);
				cleanDoubleArray(f);

				for (int l=0;l<shape_e.length;l++){
					for (int n=0;n<coordinatesreal.length;n++){
						coordinatesreal[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
						f[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(inFieldForce).getDofValues(n);//inFieldForce=2

						tksi[n]+=shape_e_der[l][0]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
						teta[n]+=shape_e_der[l][1]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);

					} //for n, coordinatesreal.length
				} // for l, shape_e.length



				normale = GeometryTools.CrossProduct(tksi, teta);
				GeometryTools.Normalize(normale);

				gksiksi=GeometryTools.ScalarProduct(tksi, tksi);
				gksieta=GeometryTools.ScalarProduct(tksi, teta);
				getaeta=GeometryTools.ScalarProduct(teta, teta);

				g = gksiksi*getaeta-gksieta*gksieta;
				//wkg = wk*Math.sqrt(g);
				wkg = wk*Math.sqrt(g)/pi8;

				//target position loops, R
				for (int i=0;i<Nnodes;i++){
					x0 = xx0[i];

					y[0]=coordinatesreal[0]-x0[0];
					y[1]=coordinatesreal[1]-x0[1];
					y[2]=coordinatesreal[2]-x0[2];

					R2 = y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
					oneOverR =1./Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+epsilon*epsilon);

					oneOverR2= oneOverR*oneOverR;
					oneOverR3= oneOverR2*oneOverR;

					G00 = (R2+2*epsilon*epsilon+y[0]*y[0])*oneOverR3;
					G01 = oneOverR3*y[0]*y[1];
					G02 = oneOverR3*y[0]*y[2];
					G11 = (R2+2*epsilon*epsilon+y[1]*y[1])*oneOverR3;
					G12 = oneOverR3*y[1]*y[2];
					G22 = (R2+2*epsilon*epsilon+y[2]*y[2])*oneOverR3;

					//f[] is the physical coords values of forces
					v_n0[0]= wkg*(G00*f[0]+G01*f[1]+G02*f[2]);
					v_n0[1]= wkg*(G01*f[0]+G11*f[1]+G12*f[2]);
					v_n0[2]= wkg*(G02*f[0]+G12*f[1]+G22*f[2]);


					vitn0x[i] += v_n0[0];
					vitn0y[i] += v_n0[1];
					vitn0z[i] += v_n0[2];

				}//end nodes
			}//end quad points
		}//end elements

		for (int i=0;i<Nnodes;i++){
			n0 = m.getFields(outFieldVelocity).getNodes(i);
			n0.setDofValueAt(0, n0.getDofValueAt(0)+vitn0x[i]);
			n0.setDofValueAt(1, n0.getDofValueAt(1)+vitn0y[i]);
			n0.setDofValueAt(2, n0.getDofValueAt(2)+vitn0z[i]);
		}
				
	}
	
	public static void computeSingleLayerVelocity_CPUReg(final Mesh m, final int inFieldForce, final int outFieldVelocity) {
		int nbProcs = m.sc.nbThreads;
		double tt = System.nanoTime();
		Thread[] productThreads = new Thread[nbProcs];
		//Dividing the domain into equals domains, except last one which contains the rest
		final int [][] bounds = new int[nbProcs][2];
		int N = m.getElements().length / nbProcs ;
		for ( int i = 0 ; i < nbProcs ; i++ )
		{
			bounds[i][0] = i * N ;
			bounds[i][1] = ( i + 1 ) * N ;
		}
		bounds[nbProcs - 1][1] = m.getElements().length;
		m.sc.messenger.message("Number of threads "+nbProcs,3);

		final int Nnodes = m.getFields().get(inFieldForce).getNodes().length;
		final double[][] xx0 = new double[m.getFields().get(inFieldForce).getNodes().length][3];
		final double[][] nvecf0 = new double[m.getFields().get(inFieldForce).getNodes().length][3];
		final double[] nscalf0 = new double[m.getFields().get(inFieldForce).getNodes().length];


		for (int i=0;i<m.getFields().get(inFieldForce).getNodes().length;i++){
			//double[] v_n0 = new double[3];
			//loop on nodes
			Elements e = m.getFields().get(inFieldForce).getNodes()[i].getVertex().getElements().get(1);
			Node n0 = m.getFields().get(inFieldForce).getNodes()[i];
			Vertex v0 = n0.getVertex();

			double[] shape;
			double[][] shape_der;
			double[][] shape_der_near;

			if (e.getVertices()[0].equals(v0)){
				shape = e.computeShapeFunctionAt(new double[]{0.0,0.0});
				shape_der=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.0,0.0});
				shape_der_near=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.1,0.1});

			}else{
				if (e.getVertices()[1].equals(v0)){
					shape = e.computeShapeFunctionAt(new double[]{1.0,0.0});
					shape_der=e.computeShapeFunctionFirstDerivativesAt(new double[]{1.0,0.0});
					shape_der_near=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.9,0.05});
				}else{
					shape = e.computeShapeFunctionAt(new double[]{0.0,1.0});
					shape_der=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.0,1.0});
					shape_der_near=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.05,0.9});
				}
			}
			//real position of vertex
			double[] x0 = new double[3];
			double[] t0 = new double[3];
			double[] t1 = new double[3];
			double[] t0_near = new double[3];
			double[] t1_near = new double[3];
			double[] f0 = new double[3];


			for (int kk=0;kk<shape.length;kk++){
				for (int jj=0;jj<x0.length;jj++){
					x0[jj]+= shape[kk]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					f0[jj]+= shape[kk]*e.getOneRing()[kk].getNodes().get(inFieldForce).getDofValues()[jj];

					t0[jj]+= shape_der[kk][0]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					t1[jj]+= shape_der[kk][1]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					t0_near[jj]+= shape_der_near[kk][0]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					t1_near[jj]+= shape_der_near[kk][1]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
				}
			}
			double[] normalenear = GeometryTools.CrossProduct(t0_near, t1_near);
			GeometryTools.Normalize(normalenear);
			double[] normale0 = GeometryTools.CrossProduct(t0, t1);
			GeometryTools.Normalize(normale0);
			if (GeometryTools.ScalarProduct(normalenear, normale0)<0.0){
				normale0 =GeometryTools.timesLambda(normale0, -1.0);
			}

			xx0[i]=x0;
			nvecf0[i]=GeometryTools.CrossProduct(normale0, f0);
			nscalf0[i]=GeometryTools.ScalarProduct(normale0, f0);
		}
		//		m.sc.messenger.message("Time for preparation  "+(System.nanoTime()-tt)/1e9,3);


		for ( int ii = 0 ; ii < nbProcs ; ii++ )
		{
			final int Ni = ii ;
			Thread threadii = new Thread( new Runnable()
			{
				int number = Ni ;


				public void run()
				{	//double Ttot =  System.nanoTime();
					Quadrature quad  = new GaussPoints12PtsOnT1();

					double G00,G01,G02,G11,G12,G22;
					
					double epsilon = m.sc.reg_epsilon;

					double[] shape_e;
					double[][] shape_e_der;
					double[] tksi= new double[3];
					double[] teta= new double[3];
					double[] coordinatesreal = new double[3];
					double[] f = new double[3];
					double wk;
					double[] normale;
					double g;
					double wkg;
					double gksiksi;
					double gksieta;
					double getaeta;
					Node n0;
					double[] x0;
					double[] nvecf;
					double f0scaln0;
					double[] v_n0 = new double[3];
					double[] vitn0x =new double[nscalf0.length];
					double[] vitn0y =new double[nscalf0.length];
					double[] vitn0z =new double[nscalf0.length];

					double[][] xx0thread =new double[nscalf0.length][3];
					double[][] nvecf0thread =new double[nscalf0.length][3];
					double[] nscalf0thread =new double[nscalf0.length];


					for (int i=0;i<Nnodes;i++){
						xx0thread[i][0]=xx0[i][0];
						xx0thread[i][1]=xx0[i][1];
						xx0thread[i][2]=xx0[i][2];

						nvecf0thread[i][0] = nvecf0[i][0];
						nvecf0thread[i][1] = nvecf0[i][1];
						nvecf0thread[i][2] = nvecf0[i][2];
						nscalf0thread[i] = nscalf0[i];
					}



					int nbquad = quad.getNumberOfQuadraturePoints();
					double[] y = new double[3];
					double pi8=8.0*Math.PI;
					double oneOverR;
					double oneOverR2;
					double oneOverR3;
					double R2;

					//for (int j=0;j<m.getElements().length;j++){
					for ( int j = bounds[number][0] ; j < bounds[number][1] ; j++ ){

						for (int k=0;k<nbquad;k++){

							wk = quad.getWeightPoint(k);
							shape_e = m.getElements(j).computeShapeFunctionAt(quad.getQuadraturePoint(k));
							shape_e_der = m.getElements(j).computeShapeFunctionFirstDerivativesAt(quad.getQuadraturePoint(k));

							cleanDoubleArray(tksi);
							cleanDoubleArray(teta);
							cleanDoubleArray(coordinatesreal);
							cleanDoubleArray(f);


							for (int l=0;l<shape_e.length;l++){
								for (int n=0;n<coordinatesreal.length;n++){
									coordinatesreal[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
									f[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(inFieldForce).getDofValues(n);

									tksi[n]+=shape_e_der[l][0]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
									teta[n]+=shape_e_der[l][1]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);

								}
							}



							normale = GeometryTools.CrossProduct(tksi, teta);
							GeometryTools.Normalize(normale);

							gksiksi=GeometryTools.ScalarProduct(tksi, tksi);
							gksieta=GeometryTools.ScalarProduct(tksi, teta);
							getaeta=GeometryTools.ScalarProduct(teta, teta);

							g = gksiksi*getaeta-gksieta*gksieta;
							//wkg = wk*Math.sqrt(g);
							wkg = wk*Math.sqrt(g)/pi8;

							///////////////
							//f = normale;
							////////////////

							double tt0 = System.nanoTime();
							for (int i=0;i<Nnodes;i++){
								x0 = xx0[i];

								//computeGlnAndRlnAt(coordinatesreal, x0, Gln, Rln);
								y[0]=coordinatesreal[0]-x0[0];
								y[1]=coordinatesreal[1]-x0[1];
								y[2]=coordinatesreal[2]-x0[2];

								R2 = y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
								oneOverR =1./Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+epsilon*epsilon);

								oneOverR2= oneOverR*oneOverR;
								oneOverR3= oneOverR2*oneOverR;
								//oneOverR3 /= (4.0*pi);
								//oneOverR /= (8.0*pi);

								G00 = (R2+2*epsilon*epsilon+y[0]*y[0])*oneOverR3;
								G01 = oneOverR3*y[0]*y[1];
								G02 = oneOverR3*y[0]*y[2];
								G11 = (R2+2*epsilon*epsilon+y[1]*y[1])*oneOverR3;
								G12 = oneOverR3*y[1]*y[2];
								G22 = (R2+2*epsilon*epsilon+y[2]*y[2])*oneOverR3;

								v_n0[0]= wkg*(G00*f[0]+G01*f[1]+G02*f[2]);
								v_n0[1]= wkg*(G01*f[0]+G11*f[1]+G12*f[2]);
								v_n0[2]= wkg*(G02*f[0]+G12*f[1]+G22*f[2]);


								//double[] vitn0x =new double[nscalf0.length];
								vitn0x[i] += v_n0[0];
								vitn0y[i] += v_n0[1];
								vitn0z[i] += v_n0[2];


							}//end nodes
							//tta+=(System.nanoTime()-tt0)/1e9;
						}//end quad points
					}//end elements

					for (int i=0;i<Nnodes;i++){

						n0 = m.getFields(outFieldVelocity).getNodes(i);
						synchronized (n0) {
							n0.setDofValueAt(0, n0.getDofValueAt(0)+vitn0x[i]);
							n0.setDofValueAt(1, n0.getDofValueAt(1)+vitn0y[i]);
							n0.setDofValueAt(2, n0.getDofValueAt(2)+vitn0z[i]);
						}
					}
				
				}				

			} ) ;

			productThreads[ii] = threadii ;
		}

		for ( int i = 0 ; i < nbProcs ; i++ )
		{
			productThreads[i].start() ;
		}

		try
		{
			for ( int i = 0 ; i < nbProcs ; i++ )
			{
				productThreads[i].join() ;
			}
		}
		catch ( InterruptedException e )
		{
			// Restore the interrupted status
			Thread.currentThread().interrupt();

		};

		m.sc.messenger.message("Time to compute velocity field on CPU "+(System.nanoTime()-tt)/1e9,3);

	}
	
	public static void computeSingleLayerVelocity_CPUCurrent(final Mesh m, final int inFieldForce, final int outFieldVelocity) {
		int nbProcs = m.sc.nbThreads;
		double tt = System.nanoTime();
		Thread[] productThreads = new Thread[nbProcs];
		//Dividing the domain into equals domains, except last one which contains the rest
		final int [][] bounds = new int[nbProcs][2];
		int N = m.getElements().length / nbProcs ;
		for ( int i = 0 ; i < nbProcs ; i++ )
		{
			bounds[i][0] = i * N ;
			bounds[i][1] = ( i + 1 ) * N ;
		}
		bounds[nbProcs - 1][1] = m.getElements().length;
		m.sc.messenger.message("Number of threads "+nbProcs,3);

		final int Nnodes = m.getFields().get(outFieldVelocity).getNodes().length;
		final double[][] xx0 = new double[m.getFields().get(outFieldVelocity).getNodes().length][3];
		final double[][] nvecf0 = new double[m.getFields().get(outFieldVelocity).getNodes().length][3];
		final double[] nscalf0 = new double[m.getFields().get(outFieldVelocity).getNodes().length];


		for (int i=0;i<m.getFields().get(outFieldVelocity).getNodes().length;i++){
			//double[] v_n0 = new double[3];
			//loop on nodes
			Elements e = m.getFields().get(outFieldVelocity).getNodes()[i].getVertex().getElements().get(1);
			Node n0 = m.getFields().get(outFieldVelocity).getNodes()[i];
			Vertex v0 = n0.getVertex();

			double[] shape;
			double[][] shape_der;
			double[][] shape_der_near;

			if (e.getVertices()[0].equals(v0)){
				shape = e.computeShapeFunctionAt(new double[]{0.0,0.0});
				shape_der=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.0,0.0});
				shape_der_near=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.1,0.1});

			}else{
				if (e.getVertices()[1].equals(v0)){
					shape = e.computeShapeFunctionAt(new double[]{1.0,0.0});
					shape_der=e.computeShapeFunctionFirstDerivativesAt(new double[]{1.0,0.0});
					shape_der_near=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.9,0.05});
				}else{
					shape = e.computeShapeFunctionAt(new double[]{0.0,1.0});
					shape_der=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.0,1.0});
					shape_der_near=e.computeShapeFunctionFirstDerivativesAt(new double[]{0.05,0.9});
				}
			}
			//real position of vertex
			double[] x0 = new double[3];
			double[] t0 = new double[3];
			double[] t1 = new double[3];
			double[] t0_near = new double[3];
			double[] t1_near = new double[3];
			double[] f0 = new double[3];


			for (int kk=0;kk<shape.length;kk++){
				for (int jj=0;jj<x0.length;jj++){
					x0[jj]+= shape[kk]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					f0[jj]+= shape[kk]*e.getOneRing()[kk].getNodes().get(inFieldForce).getDofValues()[jj];

					t0[jj]+= shape_der[kk][0]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					t1[jj]+= shape_der[kk][1]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					t0_near[jj]+= shape_der_near[kk][0]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
					t1_near[jj]+= shape_der_near[kk][1]*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];
				}
			}
			double[] normalenear = GeometryTools.CrossProduct(t0_near, t1_near);
			GeometryTools.Normalize(normalenear);
			double[] normale0 = GeometryTools.CrossProduct(t0, t1);
			GeometryTools.Normalize(normale0);
			if (GeometryTools.ScalarProduct(normalenear, normale0)<0.0){
				normale0 =GeometryTools.timesLambda(normale0, -1.0);
			}

			xx0[i]=x0;
			nvecf0[i]=GeometryTools.CrossProduct(normale0, f0);
			nscalf0[i]=GeometryTools.ScalarProduct(normale0, f0);
		}
		//		m.sc.messenger.message("Time for preparation  "+(System.nanoTime()-tt)/1e9,3);


		for ( int ii = 0 ; ii < nbProcs ; ii++ )
		{
			final int Ni = ii ;
			Thread threadii = new Thread( new Runnable()
			{
				int number = Ni ;


				public void run()
				{	//double Ttot =  System.nanoTime();
					Quadrature quad  = new GaussPoints12PtsOnT1();

					double G00,G01,G02,G11,G12,G22;
					double R00,R01,R02,R11,R12,R22;

					double[] shape_e;
					double[][] shape_e_der;
					double[] tksi= new double[3];
					double[] teta= new double[3];
					double[] coordinatesreal = new double[3];
					double[] f = new double[3];
					double wk;
					double[] normale;
					double g;
					double wkg;
					double gksiksi;
					double gksieta;
					double getaeta;
					Node n0;
					double[] x0;
					double[] nvecf;
					double[] ftilde=new double[3];
					double f0scaln0;
					double[] Rnorm=new double[3];
					double[] v_n0 = new double[3];
					double[] vitn0x =new double[nscalf0.length];
					double[] vitn0y =new double[nscalf0.length];
					double[] vitn0z =new double[nscalf0.length];

					double[][] xx0thread =new double[nscalf0.length][3];
					double[][] nvecf0thread =new double[nscalf0.length][3];
					double[] nscalf0thread =new double[nscalf0.length];


					for (int i=0;i<Nnodes;i++){
						xx0thread[i][0]=xx0[i][0];
						xx0thread[i][1]=xx0[i][1];
						xx0thread[i][2]=xx0[i][2];

						nvecf0thread[i][0] = nvecf0[i][0];
						nvecf0thread[i][1] = nvecf0[i][1];
						nvecf0thread[i][2] = nvecf0[i][2];
						nscalf0thread[i] = nscalf0[i];
					}



					int nbquad = quad.getNumberOfQuadraturePoints();
					double[] y = new double[3];
					double pi=Math.PI;
					double oneOverR;
					double oneOverR2;
					double oneOverR3;

					//for (int j=0;j<m.getElements().length;j++){
					for ( int j = bounds[number][0] ; j < bounds[number][1] ; j++ ){

						for (int k=0;k<nbquad;k++){

							wk = quad.getWeightPoint(k);
							shape_e = m.getElements(j).computeShapeFunctionAt(quad.getQuadraturePoint(k));
							shape_e_der = m.getElements(j).computeShapeFunctionFirstDerivativesAt(quad.getQuadraturePoint(k));

							cleanDoubleArray(tksi);
							cleanDoubleArray(teta);
							cleanDoubleArray(coordinatesreal);
							cleanDoubleArray(f);


							for (int l=0;l<shape_e.length;l++){
								for (int n=0;n<coordinatesreal.length;n++){
									coordinatesreal[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
									f[n]+=shape_e[l]*m.getElements(j).getOneRing(l).getNodes(inFieldForce).getDofValues(n);

									tksi[n]+=shape_e_der[l][0]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);
									teta[n]+=shape_e_der[l][1]*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);

								}
							}



							normale = GeometryTools.CrossProduct(tksi, teta);
							GeometryTools.Normalize(normale);

							gksiksi=GeometryTools.ScalarProduct(tksi, tksi);
							gksieta=GeometryTools.ScalarProduct(tksi, teta);
							getaeta=GeometryTools.ScalarProduct(teta, teta);

							g = gksiksi*getaeta-gksieta*gksieta;
							//wkg = wk*Math.sqrt(g);
							wkg = wk*Math.sqrt(g);

							///////////////
							//f = normale;
							////////////////

							double tt0 = System.nanoTime();
							for (int i=0;i<Nnodes;i++){
								x0 = xx0[i];

								//computeGlnAndRlnAt(coordinatesreal, x0, Gln, Rln);
								y[0]=coordinatesreal[0]-x0[0];
								y[1]=coordinatesreal[1]-x0[1];
								y[2]=coordinatesreal[2]-x0[2];
								
								oneOverR =1./Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);

								oneOverR2= oneOverR*oneOverR;
								oneOverR3= oneOverR2*oneOverR;
								oneOverR3 /= (4.0*pi);
								oneOverR /= (8.0*pi);

								G00 = oneOverR+y[0]*y[0]*oneOverR3*0.5;
								G01 = oneOverR3*y[0]*y[1]*0.5;
								G02 = oneOverR3*y[0]*y[2]*0.5;
								G11 = oneOverR+y[1]*y[1]*oneOverR3*0.5;
								G12 = oneOverR3*y[1]*y[2]*0.5;
								G22 = oneOverR+y[2]*y[2]*oneOverR3*0.5;

								R00 = y[0]*y[0]*oneOverR3;
								R01 = y[0]*y[1]*oneOverR3;
								R02 = y[0]*y[2]*oneOverR3;
								R11 = y[1]*y[1]*oneOverR3;
								R12 = y[1]*y[2]*oneOverR3;
								R22 = y[2]*y[2]*oneOverR3;



								nvecf = nvecf0[i];
								f0scaln0 = nscalf0[i];

								ftilde[0]=f[0]-normale[0]*f0scaln0+normale[1]*nvecf[2]-normale[2]*nvecf[1];
								ftilde[1]=f[1]-normale[1]*f0scaln0+normale[2]*nvecf[0]-normale[0]*nvecf[2];
								ftilde[2]=f[2]-normale[2]*f0scaln0+normale[0]*nvecf[1]-normale[1]*nvecf[0];
								//
								Rnorm[0]= R00*normale[0]+R01*normale[1]+R02*normale[2];
								Rnorm[1]= R01*normale[0]+R11*normale[1]+R12*normale[2];
								Rnorm[2]= R02*normale[0]+R12*normale[1]+R22*normale[2];

								v_n0[0]= wkg*(G00*ftilde[0]+G01*ftilde[1]+G02*ftilde[2]+Rnorm[1]*nvecf[2]-Rnorm[2]*nvecf[1]);
								v_n0[1]= wkg*(G01*ftilde[0]+G11*ftilde[1]+G12*ftilde[2]+Rnorm[2]*nvecf[0]-Rnorm[0]*nvecf[2]);
								v_n0[2]= wkg*(G02*ftilde[0]+G12*ftilde[1]+G22*ftilde[2]+Rnorm[0]*nvecf[1]-Rnorm[1]*nvecf[0]);


								//double[] vitn0x =new double[nscalf0.length];
								vitn0x[i] += v_n0[0];
								vitn0y[i] += v_n0[1];
								vitn0z[i] += v_n0[2];


							}//end nodes
							//tta+=(System.nanoTime()-tt0)/1e9;
						}//end quad points
					}//end elements

					for (int i=0;i<Nnodes;i++){

						n0 = m.getFields(outFieldVelocity).getNodes(i);
						synchronized (n0) {
							n0.setDofValueAt(0, n0.getDofValueAt(0)+vitn0x[i]);
							n0.setDofValueAt(1, n0.getDofValueAt(1)+vitn0y[i]);
							n0.setDofValueAt(2, n0.getDofValueAt(2)+vitn0z[i]);
						}
					}
				
				}				

			} ) ;

			productThreads[ii] = threadii ;
		}

		for ( int i = 0 ; i < nbProcs ; i++ )
		{
			productThreads[i].start() ;
		}

		try
		{
			for ( int i = 0 ; i < nbProcs ; i++ )
			{
				productThreads[i].join() ;
			}
		}
		catch ( InterruptedException e )
		{
			// Restore the interrupted status
			Thread.currentThread().interrupt();

		};

		m.sc.messenger.message("Time to compute velocity field on CPU "+(System.nanoTime()-tt)/1e9,3);

	}
	
	//utility function : puts all value of an array of double at 0.0
	private static void cleanDoubleArray(double[] arr){
		for (int i=0;i<arr.length;i++){
			arr[i]=0.0;
		}
	}
	
	public static void applyPrefixForce(Mesh m, final int inFieldForceLim, final int outFieldForce) {
		if (m.sc.nbThreads==1) {
			applyPrefixForce_1T(m, inFieldForceLim, outFieldForce);
		} else {
			System.out.println("Multi apply force has not realized yet!");
		}
	}
	
	//inFieldForceLim=5, outFieldForce=2
	public static void applyPrefixForce_1T(Mesh m, final int inFieldForceLim, final int outFieldForce) {
		double[] pos = new double[m.getFields().get(0).getNodes().length*3];
		double[] flim = new double[m.getFields().get(0).getNodes().length*3];
		double[] fv = new double[m.getFields().get(0).getNodes().length*3];
		for (int j=0;j<3;j++){
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++){
				pos[3*i+j] = m.getFields().get(4).getNodes()[i].getDofValueAt(j);
			}
		}
		//set the applied force at physical space
		for (int i=0;i<m.getFields().get(0).getNodes().length;i++) {
			flim[3*i+0] = pos[3*i+1]*pos[3*i+2];
			flim[3*i+1] = pos[3*i+2]*pos[3*i+0];
			flim[3*i+2] = pos[3*i+0]*pos[3*i+1];
		}
		m.resetField(2);
		m.resetField(5);
		for (int j=0;j<3;j++){
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++){
				m.getFields().get(inFieldForceLim).getNodes()[i].setDofValueAt(j, flim[3*i+j]);
			}
		}
		//convert to nodal values
		fv = m.convertValueAtNodeInNodalValues(flim);
		for (int j=0;j<3;j++){
			for (int i=0;i<m.getFields().get(0).getNodes().length;i++){
				m.getFields().get(outFieldForce).getNodes()[i].setDofValueAt(j, fv[3*i+j]);
			}
		}
		
	}
	
	public static Mesh prepare(){
		sc.messenger = new Messenger(sc.intparams.getParameter("Verbosity level"));
		sc.repertoire_calcul = sc.stringparams.getParameter("Computation directory");
		sc.nbThreads = sc.intparams.getParameter("Number of threads");
		
		while(sc.arch==null){
			sc.arch=sc.stringparams.getParameter("Computation of velocity field");
			if (sc.arch.equals("GPU")){
				System.out.println("Not implented yet.");
			}
		}
		
		while(sc.sparsestorage==null){
			String sparse = sc.stringparams.getParameter("Sparse matrix solver");
			if (sparse.equals("UMFPACK")){
				sc.sparsestorage="column";
			}else{
				if (sparse.equals("PARDISO")){
					sc.sparsestorage="row";
				}else{
					sc.stringparams.removeParameter("Sparse matrix solver");
					System.out.println("Typo or sparse matrix solver not implemented. Please type again.");
					System.out.println("Possible choices are : UMFPACK, PARDISO (on mesocenter) ");
				}
			}
		}
		
		while(sc.fullstorage==null){
			String full = sc.stringparams.getParameter("Full matrix solver");
			if (full.equals("UJMP")){
				sc.fullstorage="ujmp";
			}else{
				if (full.equals("LAPACK")){
					sc.fullstorage="double";
				}else{
					sc.stringparams.removeParameter("Full matrix solver");
					System.out.println("Typo or full matrix solver not implemented. Please type again.");
					System.out.println("Possible choices are : UJMP, LAPACK (on mesocenter) ");
				}
			}
		}
		
		sc.regularized = sc.intparams.getParameter("Regularized force");
		if (sc.regularized!=0) {
			sc.reg_epsilon = sc.doubleparams.getParameter("Regularized strength");
		}
		
		System.out.println("Creating mesh....");
		sc.Nsub = sc.intparams.getParameter("Number of icosahedron refinement");
		Mesh m1 = Mesh.createSphereMeshWithInterpolation(sc.Nsub);
		
		SimpleDateFormat formater = null;
		Date aujourdhui = new Date();
		formater = new SimpleDateFormat("dd-MM-yy-HH:mm:ss");
		sc.repertoire_simu = sc.repertoire_calcul.concat("SL_"+formater.format(aujourdhui));
		Field f0 = new VectorField(m1, "pos_ref");
		Field f1 = new VectorField(m1, "pos");
		Field f2 = new VectorField(m1, "force");
		Field f3 = new VectorField(m1, "velocity");
		Field f4 = new VectorField(m1, "pos_lim");
		Field f5 = new VectorField(m1, "force_lim");
		Field f6 = new VectorField(m1, "velocity_lim");
		
		for (int i = 0; i < f0.getNodes().length; i++) {
			Vertex v = f0.getNodes()[i].getVertex();
			f0.getNodes()[i].setDofValues(v.getCoordinates());
			f1.getNodes()[i].setDofValues(v.getCoordinates());
		}
		//compute the pos_lim
		m1.resetField(4);
        m1.computeLimitValueOfField1AndStoreInField2(1, 4);
        
		f2.turnActive();
		f3.turnActive();
		
		if (sc.regularized!=0) {
			sc.repertoire_simu=sc.repertoire_simu.concat("_Reg_"+String.valueOf(sc.reg_epsilon));
		}
		sc.repertoire_simu=sc.repertoire_simu.concat("/");
		//Creates the simulation directory (needs the method softobject.prepare() to have been run).
		File fb = new File(sc.repertoire_simu); 
		fb.mkdirs();
		
		//associate the mesh with the simulation context
		m1.setSimulationContext(sc);
		m1.computeCollocationMatrixNEW();
		
		System.out.println("Nb of elements "+m1.getElements().length);
		System.out.println("Nb of nodes "+m1.getFields().get(0).getNodes().length);
		
		return m1;
	}// end for prepare() method

}
