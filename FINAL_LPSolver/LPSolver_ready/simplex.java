import java.io.*;
import java.lang.*;
import java.util.StringTokenizer;

public class simplex
{
    int simpSolveFlg=0;
    addVar ad;
    public static double sTable[][] = new double[15][35];
    public static double tempMat[][] = new double[15][35];
    public static int XB[]  = new int[15];
    public static double Cobj[]  = new double[35];
    public static double RHS[]  = new double[15];
    public static double dmyRHS[]  = new double[15];
    public static double PbarJ[];
    public static double CB[] = new double[15];
    public static double Y[];
    public static double theta[];
    public static double B[][]; 
    public static double BInverse[][];
    public static double Cj_Zj[] = new double[35];
    public static int rest[] = new int[15];
    String vars[] = new String[15];
    static int slk[] = new int[35];  
    static int art[] = new int[15];  
    public static String data = new String();
    public static String errdata = new String();
    
    String tmp,tk;
    static int i,j,artCnt=0,TnoV=0,noV=0,noE=0,lineNo=0,errorFlg=0,objFlg=0,minFlg=0;

    public String cnstrNBV()
    {
	data = "\n* Range of Values for Constraint Coefficients of NBVs:-\n";
	int basFlg,brkFlg,k;
	double tmpD=0,coeff=0;
	for(k=0;k<noE;k++)
	    {
		for(i=0;i<TnoV;i++)
		    {
			basFlg = 0;
			for(j=0;j<noE;j++)
			    {
				if(XB[j] == i)
				    {
					basFlg = 1;
					break;
				    }
			    }
			if(basFlg == 1)
			    continue;
			/*now Xi is NBV */
			tmpD = 0;
			brkFlg = 0;
			for(j=0;j<noE;j++)
			    {
				if(j == k)
				    {
					coeff = Y[j];
					if(coeff == 0)
					    {
						brkFlg = 1;
						tmpD = 0;
						break;
					    }
				    }
				else
				    tmpD += Y[j] * sTable[j][i];
			    }
			tmpD = Cobj[i] - tmpD;
			if(brkFlg == 0)
			    {
				if(coeff < 0)
				    data += "\ta"+(k+1)+(i+1)+" <= "+((-1)*tmpD/(-1 *coeff))+"\n";
				else
				    data += "\ta"+(k+1)+(i+1)+" >= "+(tmpD/coeff)+"\n";
			    }
		    }
		data += "\n";
	    }
	return data;
    }
    
    public String cnstrRHS()
    {
	int basFlg,k;
	double tmpD=0,coeff=0,lb,ub;
	data = "\n* Range of Values for RHS values of Constraints:-\n";
	for(j=0;j<noE;j++)
	    {
		lb=-999999;
		ub=999999;
		for(k=0;k<noE;k++)
		    {
			tmpD =0;
			basFlg = 0;
			for(i=0;i<noE;i++)
			    {
				if(i == j)
				    {
					coeff = BInverse[k][i];
					if(coeff == 0)
					    {
						basFlg = 1;
						tmpD = 0;
						break;
					    }
				    }
				else
				    tmpD += BInverse[k][i] * dmyRHS[i];
			    }
			if(basFlg == 0)
			    {
				if(coeff < 0 )
				    {
					data += "\tb"+(j+1)+" <= "+(tmpD/(-1 *coeff))+"\n";
					if((tmpD/(-1 *coeff)) < ub)
					    ub = (tmpD/(-1 *coeff));
				    }
				else
				    {
					data += "\tb"+(j+1)+" >= "+((-1)*tmpD/coeff)+"\n";
					if(((-1)*tmpD/coeff) > lb)
					    lb = ((-1)*tmpD/coeff);
				    }
			    }
		    }
		data += "\ti.e. "+lb+" <= b"+(j+1)+" <= "+ub+"\n";
	    }
	
	return data;
    }

    public String objBV()
    {
	int brkFlg,basFlg,k,l;
	double tmpD=0,coeff=0,lb,ub;
	data = "\n* Range of Values for Objective Function Coefficients of BVs:-\n";
	for(i=0;i<noE;i++)
	    {
		lb=-99999999;
		ub=99999999;
		for(j=0;j<TnoV;j++)
		    {
			basFlg = 0;
			for(l=0;l<noE;l++)
			    if(XB[l] == j)
				{
				    basFlg = 1;
				    break;
				}
			if(basFlg == 1)
			    {
				basFlg = 0;
				continue;
			    }
			/*now Xj is NBV */
			calcPbarJ(j);
			tmpD = 0;
			brkFlg = 0;
			for(k=0;k<noE;k++)
			    {
				if(k == i)
				    {
					coeff = PbarJ[k];
					if(coeff == 0)
					    {
						brkFlg = 1;
						tmpD = 0;
						break;
					    }
				    }
				else
				    tmpD += CB[k] * PbarJ[k];
			    }
			tmpD = Cobj[j] - tmpD;
			if(brkFlg == 0)
			    {
				if(coeff < 0)
				    {
					data += "\tC"+(XB[i]+1)+" <= "+((-1)*tmpD/(-1 *coeff))+"\n";
					if(((-1)*tmpD/(-1 *coeff)) < ub)
					    ub = ((-1)*tmpD/(-1 *coeff));
				    }
				else
				    {
					data += "\tC"+(XB[i]+1)+" >= "+(tmpD/coeff)+"\n";
					if((tmpD/coeff) > lb)
					    lb = (tmpD/coeff);
				    }
			    }
		    }
		data += "\ti.e. "+lb+" <= C"+(XB[i]+1)+" <= "+ub+"\n";
	    }
	return data;
    }

    public String objNBV()
    {
	int basFlg,k;
	double tmpD=0;
	data = "\n* Range of Values for Objective Function Coefficients of NBVs:-\n";
	for(i=0;i<TnoV;i++)
	    {
		basFlg = 0;
		for(j=0;j<noE;j++)
		    {
			if(XB[j] == i)
			    {
				basFlg = 1;
				break;
			    }
		    }
		if(basFlg == 1)
		    continue;
		/*now Xi is NBV */
		calcPbarJ(i);
		tmpD = 0;
		for(k=0;k<noE;k++)
		    {
			tmpD += CB[k] * PbarJ[k];
			//System.out.println("XB:"+CB[k]);
		    }
		data += "\tC"+(i+1)+" <= "+tmpD+"\n";
	    }
	return data;
    }
    
    public String addNewVar(addVar ad,simplex s1)
    {
	double d1 = 1;
	double maxCjZj;
	int cngFlg = 0;
	this.ad = ad;
	vars[TnoV] = ad.vname;
	double d = Double.parseDouble(ad.restr);
	
	data = "";
	if(d == 0) // unrestricted
	    {
		for(j=0;j<2;j++)
		    {
			vars[TnoV] = ad.vname + (j+1);
			Cobj[TnoV] = d1 * Double.parseDouble(ad.objf);
			tk="";
			int i=0;
			double Znew = 0;
			StringTokenizer sti = new StringTokenizer(ad.cnstr,",");
			while(sti.hasMoreTokens())
			    {
				tk = sti.nextToken();
				sTable[i][TnoV] = d1 * Double.parseDouble(tk);
				i++;
			    }
			TnoV++;
			calcPbarJ(TnoV-1);
			for(i=0;i<noE;i++)
			    Znew += CB[i]*PbarJ[i];
			if((Cobj[TnoV-1]-Znew) <= 0)
			    {
				data += "C"+TnoV+"-Z"+TnoV+" = "+(Cobj[TnoV-1]-Znew)+" <= 0 i.e. No Change in Solution\n";
			    }
			else
			    {
				cngFlg = 1;
			    }
			d1 = -1;
		    }
		if(cngFlg == 1)
		    {
			data ="";
			s1.simplexMethod(s1);
		    }
		else
		    {
			TnoV -= 2;
		    }
	    }
	else
	  {
		StringTokenizer sti = new StringTokenizer(ad.cnstr,",");
		d1 = d;
		Cobj[TnoV] = d1 * Double.parseDouble(ad.objf);
		tk="";
		int i=0;
		double Znew = 0;
		while(sti.hasMoreTokens())
		    {
			tk = sti.nextToken();
			sTable[i][TnoV] = d1 * Double.parseDouble(tk);
			i++;
		    }
		TnoV++;
		calcPbarJ(TnoV-1);
		for(i=0;i<noE;i++)
		    Znew += CB[i]*PbarJ[i];
		maxCjZj = Cobj[TnoV-1]-Znew;
		if(maxCjZj <= 0)
		    {
			data = "C"+TnoV+"-Z"+TnoV+" = "+maxCjZj+" <= 0 i.e. No Change in Solution\n";
			TnoV--;
		    }
		else
		    {
			data ="";
			s1.simplexMethod(s1);
		    }
	  }
	return data;
    }
    
    public void simplexMethod(simplex s1)
    {
	int i,j;
	theta = new double[noE]; 
	B = new double[noE][noE];
	PbarJ  = new double[noE];
	BInverse = new double[noE][noE];
	Y  = new double[noE];
	int dmyXB[]  = new int[15];
	int iteration=0,dmFlg = 0;
	for(i=0;i<noV;i++)
	    {
		//System.out.println(i+"///rest-"+rest[i]+"-"+vars[i]);
		if(rest[i] == -1) // <=0
		    {
			Cobj[i] *= -1;
			for(j=0;j<noE;j++)
			    sTable[j][i] *= -1;
		    }
		if(rest[i] == 0) //unrestricted
		    {
			Cobj[TnoV] = Cobj[i]*(-1);
			vars[TnoV] = "s"+(TnoV+1);
			for(j=0;j<noE;j++)
			    sTable[j][TnoV] = sTable[j][i]*(-1);
			TnoV++;
		    }
	    }
	System.out.print("Simplex Method:\n\nMAX Z = ");
	data+="\nMAX Z = ";
	for(i=0;i<TnoV-1;i++) // Cobj
	    {
		System.out.print("("+Cobj[i]+vars[i]+")+");
		data += "("+Cobj[i]+vars[i]+")+";
	    }
	System.out.println("("+Cobj[i]+vars[i]+")\n");
	System.out.print("Obj:");
	data += "("+Cobj[i]+vars[i]+")\n\nObj:";
	for(i=0;i<TnoV;i++) // Cobj
	    {
		System.out.print(" | "+Cobj[i]);
		data += " | "+Cobj[i];
	    }
	int alterFlg = 0,maxJ=0;
	while(true)
	    {
		iteration++;
		double min=99999999,max=0;
		int minI=-1,k,basicFlg,endFlg=0;
		for(i=0;i<noE;i++) 
		    {
			System.out.println("\t"+XB[i]);
		    }
		//System.out.println("Matrix B:");
		for(j=0;j<noE;j++) // fill the matrix B
		    {
			for(i=0;i<noE;i++)
			    {
				B[i][j] = sTable[i][XB[j]];
				//System.out.print("\t"+B[i][j]);
			    }
			//System.out.println();
		    }
		//System.out.println("Matrix CB:");
		for(i=0;i<noE;i++) // fill matrix CB
		    {
			CB[i] = Cobj[XB[i]];
	 		//System.out.println("\t"+CB[i]);
		    }
				    
		for(j=0;j<TnoV;j++) // calculate Cj-Zj
		    {
			basicFlg = 0;
			for(k=0;k<noE;k++)
			    {
				if(j == XB[k]) // if Xj is in basic
				    {
					Cj_Zj[j] = 0;
					basicFlg = 1;
					break;
				    }
			    }
			if(basicFlg == 1)
			    continue;
			s1.calcPbarJ(j);
			
			for(i=0;i<noE;i++) // fill Temp Matrix
			    tempMat[i][j] = PbarJ[i];
			double Zj;
			int m,n;
			for(m=0;m<noE;m++) // calc Y matrix
			    {
				Zj = 0;
				for(n=0;n<noE;n++)
				    {
					Zj += CB[n]*BInverse[n][m];
				    }
				Y[m] = Zj;
				
			    }
			Zj=0;
			for(i=0;i<noE;i++) //sTable[i][j] == P[j]
			    Zj += Y[i]*sTable[i][j];
			Cj_Zj[j] = Cobj[j] - Zj;
			if((Cj_Zj[j] > 0) && (Cj_Zj[j] > max))
			    {
				max = Cj_Zj[j];
				maxJ = j;
				endFlg = 1;
			    }
		    }
		/////////
		s1.calcPbarJ(maxJ);
			//min = RHS[0]/PbarJ[0];
		for(i=0;i<noE;i++) // calc Theta
		    {
			theta[i] = RHS[i]/PbarJ[i];
			if(theta[i] >= 0 && theta[i] < min)
			    {
				min = theta[i];
				minI = i;
			    }
		    }
		for(i=0;i<noE;i++) // fill Temp matrix
		    {
			j = XB[i];
			for(k=0;k<noE;k++)
			    {
				if(XB[k] == j)
				    tempMat[k][j] = 1;
				else
				    tempMat[k][j] = 0;
			    }
		    }
			
		System.out.println("\nSimplex Table:"+iteration);
		data += "\n---------------------------------------------------------------------------------------------------------------------------------";
		data += "\nSimplex Table step: "+iteration;
		data += "\nCB | XB";
		for(i = 0; i < TnoV; i++)
		    {
			data += "  | "+vars[i];
		    }
		    data += " | RHS | theta";
		    data += "\n---------------------------------------------------------------------------------------------------------------------------------\n";
		
		for(j=0;j<noE;j++)
		    {
			System.out.print(CB[j]);
			System.out.print(" | "+vars[XB[j]]);
			data += CB[j]+" | "+vars[XB[j]];
			for(i=0;i<TnoV;i++)
			    {
				System.out.print(" | "+tempMat[j][i]);
				data += " | "+tempMat[j][i];
			    }
			System.out.print(" | "+RHS[j]);
			data += " | "+RHS[j];
			if(endFlg == 1)
			    {
				System.out.print(" | "+theta[j]);
				data += " | "+theta[j];
			    }
			data +="\n";
			System.out.println();
		    }
		System.out.println("----------------------------------------------------------");
		System.out.print("Cj_Zj:  ");
		data += "---------------------------------------------------------------------------------------------------------------------------------";
		data += "\nCj_Zj:  ";
		
		for(i=0;i<TnoV;i++) //print Cj_Zj
		    {
			System.out.print(" | "+Cj_Zj[i]);
			data += " | "+Cj_Zj[i];
		    }
		double Zj=0;
		
		for(i=0;i<noE;i++)
		    Zj += CB[i]*RHS[i];
		System.out.println(" | "+Zj);
		data += " | "+Zj;
		
		if(minI == -1) //All theta < 0 
		    {
			System.out.print("\n\n* As Pbar("+maxJ+") contains all elements <=0, Solution Space is UNBOUNDED.\n");
			data +=  "\n\n* As PbarJ of "+vars[maxJ]+" contains all elements <=0, Solution Space is UNBOUNDED.\n";
			if(Cobj[maxJ] < 0)
			    {
				data +=  "\n* As Objective Coeff. of "+maxJ+" is -ve, Objective Value also UNBOUNDED\n";
				System.out.print("\n* As Objective Coeff. of "+maxJ+" is -ve, Objective Value also UNBOUNDED\n");
			    }
			return;
		    }
		
		
		
		//System.out.println("XB change:"+XB[minI]+" at "+minI);
		if(endFlg == 0)
		    {
			data += "\n\n***************************************\nSolution:  ";
			if(minFlg == 1)
			    {
				System.out.println("\nZ = "+(-1)*Zj);
				data += "Z = "+(-1)*Zj;
			    }
			else
			    {
				System.out.println("\nZ = "+Zj);
				data += "Z = "+Zj;
			    }
			data += "\n***************************************\n";
			System.out.println(artCnt);
			for(i=0;i<noE;i++)
			    {
				for(j=0;j<artCnt;j++)
				    {
					if(art[j] == XB[i])
					    {
						System.out.println(vars[XB[i]]);
						data += "\n* As an artificial variable "+vars[XB[i]]+" is present, Solution is PRIMAL INFEASIBLE.\n";
						System.out.println("\n* As an artificial variable "+vars[XB[i]]+" is present, Solution is PRIMAL INFEASIBLE.\n");
						break;
					    }
				    }
			    }
			int zFlg = 0;
			alterFlg = 0;
			for(i=0;i<TnoV;i++) // Alternate Solution
			    {
				zFlg = 0;
				if(Cj_Zj[i] == 0)
				    {
					for(j=0;j<noE;j++)
					    {
						if(XB[j] == i)
						    {
							zFlg = 1;
							break;
						    }
					    }
					if(zFlg == 1)
					    continue;
					else
					    {
						if(dmFlg == 0)
						    {
							for(j=0;j<noE;j++)
							    dmyXB[j] = XB[j];
							dmFlg = 1;
						    }
						alterFlg = 1;
						maxJ = i;
						XB[minI] = maxJ;
						System.out.println(maxJ+" > "+minI+" var:"+vars[i]);
						s1.calcPbarJ(maxJ);
						int t;
						//min = 999999;
						for(t=0;t<noE;t++) // calc Theta
						    {
							theta[t] = RHS[t]/PbarJ[t];
							//if(theta[t] >= 0 && theta[t] < min)
							//  {
							//min = theta[t];
							//minI = t;
							//  }
						    }
						for(t=0;t<noE;t++) // set new RHS
						    {
							if(t == minI)
							    RHS[t] = theta[t];
							else
							    RHS[t] = RHS[t] - (theta[minI]*PbarJ[t]);
						    }
						
						break;
					    }
				    }
			    }
			if(alterFlg == 0)
			    break;
			else
			    {
				int flg = 0;
				for(i=0;i<noE;i++)
				    {
					//System.out.println("XB:"+XB[i]+" dmyXB:"+dmyXB[i]);
					if(XB[i]!=dmyXB[i])
					    {
						flg = 1;
						break;
					    }
				    }
				if(flg == 0)
				    break;
				data += "\n* Alternate Solution:\n";
			    }
		    }
		else
		    {
			for(i=0;i<noE;i++) // set new RHS
			    {
				if(i == minI)
				    RHS[i] = theta[i];
				else
				    RHS[i] = RHS[i] - (theta[minI]*PbarJ[i]);
			    }
			XB[minI] = maxJ; // update XB
		    }
	    }
    }
    
    public void calcPbarJ(int j)
    {
	inverseMatB();
	double tmp;
	int m,n;
	for(m=0;m<noE;m++) // calc PbarJ
	    {
		tmp = 0;
		for(n=0;n<noE;n++)
		    {
			tmp += BInverse[m][n]*sTable[n][j];
		    }
		PbarJ[m] = tmp;
		//System.out.println("P: "+PbarJ[m]);
	    }
    }
    void inverseMatB()
    {
	int st_vrs=B.length, st_stolp=B[0].length;
	double[][]old = new double[st_vrs][st_stolp*2];
	double[][]new1 = new double[st_vrs][st_stolp*2];
	int i , j ,k ;
	
	for (i=0; i < st_vrs; i++)
	    {//ones vector                                                                                                             
		for (k=0; k < st_stolp*2; k++)
		    {
			if (k-i == st_vrs)
			    old[i][k] = 1;
			if(k < st_stolp)
			    old[i][k] = B[i][k];
		    }
	    }
	for(i=0;i<st_vrs;i++)
	    {  //zeros below the diagonal  
		for(j=0;j<st_vrs;j++)
		    {
			for(k=0;k<st_stolp*2;k++)
			    {
				if (i==j)
				    new1[i][k]=old[i][k]/old[i][i];
				else
				    new1[j][k]=old[j][k];
			    }
		    }
		old=prepisi(new1);
		for(j=i+1;j<st_vrs;j++)
		    {
			for(k=0;k<st_stolp*2;k++)
			    {
				new1[j][k]=old[j][k]-old[i][k]*old[j][i];
			    }
		    }
		old=prepisi(new1);
	    }
	//zeros above the diagonal 
        for(k=st_stolp-1;k>0;k--)
	    {
		for(i=k-1;i>=0;i--)
		    {
			for(int s1=0; s1<st_stolp*2; s1++)
			    {
				new1[i][s1]=old[i][s1]-old[k][s1]*old[i][k];
			    }
		    }
		old=prepisi(new1);
	    }
	
        for(i=0;i<st_vrs;i++)
	    {//rigt part of matrix is invers                                                                                               
		for (k=st_stolp;k<st_stolp*2;k++)
		    {
			BInverse[i][k-st_stolp]=new1[i][k];
		    }
	    }
	
    }
    public double[][] prepisi(double[][]in)
    {
	double[][]out=new double[in.length][in[0].length];
	for(int v=0;v<in.length;v++)
	    {
	    for(int s=0;s<in[0].length;s++)
		{
		out[v][s]=in[v][s];
	    }
	}
	return out;
    }

    void getDual()
    {
	System.out.println("\nDUAL:");
	data = "";
	if(minFlg == 0)
	    {
		System.out.print("MIN w=");
		data += "MIN w=";
	    }
	else
	   {
	       for(i=0;i<noV;i++)
		    Cobj[i] *= -1;
	       System.out.print("MAX w=");
	       data += "MAX w=";
	   } 
	for(i=0;i<noE;i++)
	    {
		String data1="";
		if((RHS[i]>=0)&&(i != 0))
		    {
			if(RHS[i] != 1.0)
			    data1 += "+"+RHS[i];
			else
			    data1 += "+";
		    }
		else
		    {
			if(RHS[i] != 1.0)
			    data1 += ""+RHS[i];
		    } 
		data1 += "y"+(i+1);
		System.out.print(data1);
		data += data1;
	    }
	data += "\n\nsubject_to:\n";
	System.out.println("\n\nsubject_to:\n");
	for(i=0;i<noV;i++)
	    {
		System.out.print("\t");
		data += "\t";
		for(j=0;j<noE;j++)
		    {
			String data1 = "";
			if((sTable[j][i]>=0)&&(j != 0))
			    {
				if(sTable[j][i] != 1.0)
				    data1 += "+"+sTable[j][i];
				else
				    data1 += "+";
			    }
			else
			    {
				if(sTable[j][i] != 1.0)
				    data1 += ""+sTable[j][i];
			    }
			data1 += "y"+(j+1);
			System.out.print(data1);
			data += data1;
		    }
		if(rest[i] == 0)
		    {
			System.out.print(" = ");
			data += " = ";
		    }
		else if(rest[i] == -1)
		    {
			System.out.print(" <= ");
			data += " <= ";
		    }
		else if(rest[i] == 1)
		    {
			System.out.print(" >= ");
			data += " >= ";
		    }
		
		System.out.println(Cobj[i]);
		data += Cobj[i]+"\n";
	    }
	System.out.print("\nrestrictions:  ");
	data += "\nrestrictions:  ";
	for(i=0;i<noE;i++)
	    {
		if(slk[i] == 0)
		    {
			System.out.print("y"+(i+1)+"-unrestricted, ");
			data += "y"+(i+1)+"-unrestricted, ";
		    }
		else if(slk[i] == -1)
		    {
			System.out.print("y"+(i+1)+"<=0, ");
			data += "y"+(i+1)+"<=0, ";
		    }
		else if(slk[i] == 1)
		    {
			System.out.print("y"+(i+1)+">=0, ");
			data += "y"+(i+1)+">=0, ";
		    }
	    }
	System.out.println();
    }

    public String callSimplex(String fname, String sORd,simplex s1)
    {
	//simplex s1 = new simplex();
	int eqFlg=0,declFlg=0,consFlg=0,artFlg=0;
	for(i=0;i<15;i++)
	    {
		for(j=0;j<35;j++)
		    {
			sTable[i][j] = 0;
			tempMat[i][j] = 0;
			Cobj[j] = 0;
			Cj_Zj[j] = 0;
			slk[j] = 0;
		    }
		CB[i] = 0;
		XB[i] = 0;
		RHS[i] = 0;
		rest[i] = 0;
		vars[i] = "";
	    }
	try
	    {
		TnoV=0;noV=0;noE=0;lineNo=0;errorFlg=0;objFlg=0;minFlg=0;
		artCnt = 0;
		data = "";
		errdata = "ERRORS:\n";
		tmp = "";
		tk= "";
		
		FileReader fr = new FileReader(fname);
		BufferedReader br = new BufferedReader(fr);
		String line;
		
		while((line = br.readLine()) != null) 
		    {
			lineNo++;
			//System.out.println("###line: "+line);
			if(line.equals("\n"))
			    continue;
			String lower = line.toLowerCase();
			if(declFlg == 0)
			    {
				if(lower.indexOf("declarations:")!= -1)
				    {
					declFlg = 1;
					s1.parseDeclaration(line);
				    }
			    }
			else if(lower.indexOf("max:")!= -1)
			    {
				objFlg = 1;
				System.out.println(line);
				s1.parseObjective(line);
			    }
			else if(lower.indexOf("min:")!= -1)
			    {
				minFlg = 1;
				objFlg = 1;
				s1.parseObjective(line);
				for(i=0;i<noV;i++)
				    Cobj[i] *= (-1);
			    }
			else if(eqFlg == 0)
			    {
				if(lower.indexOf("subject_to:")!= -1)
				    eqFlg = 1;
			    }
			else if(lower.indexOf("restrictions:")!= -1)
			    {
				consFlg = 1;
				s1.parseConstraints(line);
			    }
			
			else if((eqFlg == 1)&&(declFlg == 1)&&(consFlg==0))
			    {
				s1.parseEquation(line);
				//if(noE == noV)
				//  break;
			    }//else
		    }//while((line = br.readLine()) != null)
		fr.close();
		if(declFlg == 0)
		    {
			System.out.println("Syntax Error- There must be Declarations: Statement");
			errdata += "Syntax Error- There must be Declarations: Statement\n";
			errorFlg = 1;
		    }
		if(objFlg == 0)
		    {
			errdata += "Syntax Error- There must be Obj Function Statement\n";
			System.out.println("Syntax Error- There must be Obj Function Statement");
			errorFlg = 1;
		    }
		if(eqFlg != 1)
		    {
			errdata += "Syntax Error- There must be subject_to: Statement\n";
			System.out.println("Syntax Error- There must be subject_to: Statement");
			errorFlg = 1;
		    }
		if(consFlg != 1)
		    {
			errdata += "Syntax Error- There must be restrictions: Statement\n";
			System.out.println("Syntax Error- There must be restrictions: Statement");
			errorFlg = 1;
		    }
		/*if(noE<noV)
		    {
			errdata += "Error- No. of Equations < No. of Variables\n";
			System.out.println("Error- No. of Equations < No. of Variables");
			errorFlg = 1;
			}*/
		if(errorFlg == 1)
		    {
			return errdata;
		    }
		int k;
		j = noV-1;
		artCnt = 0;
		for(i=0;i<noE;i++)
		    {
			j++;
			if((slk[i] == 1)||(slk[i] == 0))
			    for(k=0;k<noE;k++)
				{
				    if(k==i)
					{
					    if(slk[i] == 0)
						art[artCnt++] = j;
					    sTable[k][j] = 1;
					    XB[i] = j;
					    if(slk[i] == 0)
						Cobj[j] = -999999;
					    CB[i] = Cobj[j];
					}
				    else
					sTable[k][j] = 0;
				}
			else if(slk[i] == -1)
			    {
				for(k=0;k<noE;k++)
				    {
					if(k==i)
					    {
						sTable[k][j] = -1;
						sTable[k][j+1] = 1;
						art[artCnt++] = j+1;
						XB[i] = j+1;
						Cobj[j+1] = -999999;
						CB[i] = -999999;
						
					    }
					
					else
					    {
						sTable[k][j] = 0;
						sTable[k][j+1] = 0;
					    }
				    }
				j++;
			    }
		    }
		if(artFlg == 1)
		    {
			System.out.println("NOTE: M = 999999 used as coeff. of artificial vars.");
			data += "\nNOTE: M = 999999 used as coeff. of artificial vars.\n";
		    }
		TnoV += noV;
		//s1.printMatrix();
		if(sORd.equals("simplex"))
		    {
			for(i=0;i<noE;i++)
			    dmyRHS[i] = RHS[i];
			//j=1;
			for(i=noV;i<TnoV;i++)
			    {
				vars[i] = "s"+(i+1);
				//j++;
			    }
			s1.simplexMethod(s1);
		    }
		else
		    s1.getDual();
		
	    }
	catch(Exception e)
	    {
		System.out.println("Error Line-"+lineNo+": "+e);
		errdata += "Error Line-"+lineNo+": "+e+"\n";
		errorFlg = 1;
	    }
	return data;
    }
    
    void parseConstraints(String line)
    {
	int flg=0,noCV=0,incFlg=0,lgeqFlg=0;
	
	String tk;
	String lower = line.toLowerCase();
	String decl;
	decl = line.substring(lower.indexOf("restrictions:")+13);
	if(decl.indexOf(";")== -1)
	    {
		errorFlg = 1;
		System.out.println("In constraints, ; is must");
		errdata += "In constraints, ; is must\n";
	    }
	StringTokenizer sti = new StringTokenizer(decl,";");
	
	String tmp = "";
	while(sti.hasMoreTokens())
	    {
		tk = sti.nextToken();
		if(tk.indexOf("<=")!= -1)
		    {
			StringTokenizer stl = new StringTokenizer(tk,"<=");
			while(stl.hasMoreTokens())
			    {
				String t1 = stl.nextToken();
				StringTokenizer stle = new StringTokenizer(t1,", ");
				while(stle.hasMoreTokens())
				    {
					String t2 = stle.nextToken();
					for(i=0;i<noV;i++)
					    {
						if(t2.equals(vars[i]))
						    {
							rest[i] = -1;
							//System.out.println(i+"rest--->-1"+vars[i]);
						    }
					    }
				    }
				//0
			    }
		    }
		else if(tk.indexOf(">=")!= -1)
		    {
			StringTokenizer stl = new StringTokenizer(tk,">=");
			while(stl.hasMoreTokens())
			    {
				String t1 = stl.nextToken();
				StringTokenizer stle = new StringTokenizer(t1,", ");
				while(stle.hasMoreTokens())
				    {
					String t2 = stle.nextToken();
					for(i=0;i<noV;i++)
					    {
						if(t2.equals(vars[i]))
						    {
							rest[i] = 1;
						    }
					    }
				    }
				//0
			    }
		    }
	    }
	//for(i=0;i<noV;i++)
	//  System.out.println(rest[i]);
		
    }

    
    void parseObjective(String line)
    {
	int flg=0,noCV=0,incFlg=0,lgeqFlg=0;
	String coVar[] = new String[15];
	String tk;
	String lower = line.toLowerCase();
	String decl;
	if(minFlg == 1)
	    decl = line.substring(lower.indexOf("min:")+4);
	else
	    decl = line.substring(lower.indexOf("max:")+4);
	StringTokenizer sti = new StringTokenizer(decl,"+- ",true);
	//System.out.println(noE+1+": "+line);
	String tmp = "";
	while(true)
	    {
		flg = 0;
		if(sti.hasMoreTokens())
		    tk = sti.nextToken();
		else
		    break;
		if(tk.equals(" "))
		    continue;
		if(tk.equals("-"))
		    {
			if(tmp.equals(""))
			    tmp += tk; 
			else
			    {
				errdata += "Line-"+lineNo+": Syntax Error\n";
				System.out.println("Line-"+lineNo+": Syntax Error");
				//System.exit(0);
			    }
			if(sti.hasMoreTokens())
			    tk = sti.nextToken();
			else 
			    break;
		    }
		while((!(tk.equals("+")))&&(!(tk.equals("-"))))
		    {
			if(!(tk.equals(" ")))
			    {    
				flg = 1;
				tmp += tk;
			    }
			if(sti.hasMoreTokens())
			    tk = sti.nextToken();
			else
			    break;
		    }
		if(flg == 0)
		    {
			errdata += "Line-"+lineNo+": Syntax Error\n";
			System.out.println("Line-"+lineNo+": Syntax Error");
			errorFlg = 1;
		    }
		if(tk.equals("+"))
		    {
			coVar[noCV++] = tmp;
			//System.out.println(coVar[noCV-1]);
			tmp = "";
			continue;
		    }
		if(tk.equals("-"))
		    {
			coVar[noCV++] = tmp;
			//System.out.println(coVar[noCV-1]);
			tmp = "-";
			continue;
		    }
	    }//while(sti.hasMoreTokens())
	if(!tmp.equals(""))
	    {
		
		coVar[noCV++] = tmp;
		//System.out.println(coVar[noCV-1]);
	    }
	for(i=0;i<noCV;i++)
	    {
		for(j=0;j<noV;j++)
		    {
			if(coVar[i].endsWith(vars[j]))//Validate coeff*Var
			    {
				if(coVar[i].equals(vars[j]))
				    {
					Cobj[j] += 1;
					incFlg = 1;
				    }
				else
				    {
					tmp = coVar[i].substring(0,(coVar[i].lastIndexOf(vars[j])));
					if(tmp.equals("-"))
					    tmp = "-1";
					//System.out.println("tmp:"+tmp);
					Cobj[j] += Double.parseDouble(tmp);
					//incFlg=1;
				    }
				break;
			    }
		    }
	    }
    }

    void printMatrix() // Print AX=B
    {
	
	System.out.println("noV:"+noV);
	for(i=0;i<noE;i++)
	    System.out.println("slk: "+slk[i]);
	System.out.print("Declarations:");
	for(i=0;i<noV;i++)
	    System.out.print(" "+vars[i]);
	System.out.print("\n\nMAX Z =");
	for(i=0;i<noV;i++)
	    System.out.print(" "+Cobj[i]+vars[i]);
	System.out.println("\n\nAX = B :");
	for(i=0;i<noE;i++)
	    {
		System.out.print("[ ");
		for(j=0;j<TnoV;j++)
		    System.out.print(sTable[i][j]+" ");
		System.out.print("] [ "+vars[i]);
		System.out.println(" ] = [ "+RHS[i]+" ]");
	    }
	System.out.println("XB:");
	for(i=0;i<noE;i++)
	    System.out.println(""+XB[i]);
	System.out.println("CB:");
	for(i=0;i<noE;i++)
	    System.out.println(""+CB[i]);
    }
    void parseDeclaration(String line)
    {
	String lower = line.toLowerCase();
	String decl = line.substring(lower.indexOf("declarations:")+13);
	StringTokenizer std = new StringTokenizer(decl,", ");
	while(std.hasMoreTokens()) 
	    {
		tk = std.nextToken();
		vars[noV++] = tk;
	    }
	//System.out.println("noV:"+noV);
    }

    void parseEquation(String line)
    {
	int flg=0,noCV=0,incFlg=0,lgeqFlg=0,flipFlg=0;
	String coVar[] = new String[15];
	String tk;
	StringTokenizer sti = new StringTokenizer(line,"+-<>= ",true);
	//System.out.println(noE+1+": "+line);
	String tmp = "";
	while(sti.hasMoreTokens()) 
	    {
		flg = 0;
		tk = sti.nextToken();
		if(tk.equals(" "))
		    continue;
		if(tk.equals("-"))
		    {
			if(tmp.equals(""))
			    tmp += tk; 
			else
			    {
				System.out.println("Line-"+lineNo+": Syntax Error");
				errdata += "Line-"+lineNo+": Syntax Error\n";
				//System.exit(0);
			    }
			tk = sti.nextToken();
		    }
		while((!(tk.equals("+")))&&(!(tk.equals("-")))&&(!(tk.equals("=")))&&(!(tk.equals(">")))&&(!(tk.equals("<"))))
		    {
			if(!(tk.equals(" ")))
			    {    
				flg = 1;
				tmp += tk;
			    }
			tk = sti.nextToken();
		    }
		if(flg == 0)
		    {
			errdata += "Line-"+lineNo+": Syntax Error\n";
			System.out.println("Line-"+lineNo+": Syntax Error");
			errorFlg = 1;
		    }
		if(tk.equals("+"))
		    {
			coVar[noCV++] = tmp;
			//System.out.println(coVar[noCV-1]);
			tmp = "";
			continue;
		    }
		if(tk.equals("-"))
		    {
			coVar[noCV++] = tmp;
			//System.out.println(coVar[noCV-1]);
			tmp = "-";
			continue;
		    }
		String strLGEQ="";
		if(tk.equals("<"))
		    {
			while(sti.hasMoreTokens())
			    {
				tk =sti.nextToken();
				if(tk.equals("="))
				    {
					lgeqFlg = 1;
					strLGEQ = "<=";
					break;
				    }
				if(tk.equals(" "))
				    continue;
				else
				    {
					errdata += " <, ! <= ERROR\n";
					System.out.println(" <, ! <= ERROR");
					errorFlg = 1;
				    }
			    }
		    }
		else if(tk.equals(">"))
		    {
			while(sti.hasMoreTokens())
			    {
				tk =sti.nextToken();
				if(tk.equals("="))
				    {
					lgeqFlg = 1;
					strLGEQ = ">=";
					break;
				    }
				if(tk.equals(" "))
				    continue;
				else
				    {
					errdata += " >, ! >= ERROR\n";
					System.out.println(">, ! >= ERROR");
					errorFlg = 1;
				    }
			    }
		    }
		else if(tk.equals("="))
		    {
			lgeqFlg = 1;
			strLGEQ = "=";
		    }
		if(lgeqFlg == 1)
		    {
			coVar[noCV++] = tmp;
			//System.out.println(coVar[noCV-1]);
			coVar[noCV] = strLGEQ;
			tmp = "";
			while(sti.hasMoreTokens())
			    {
				tk = sti.nextToken();
				if(!tk.equals(" "))
				    tmp += tk;
			    }
			RHS[noE] = Double.parseDouble(tmp);
			if(RHS[noE] < 0)
			    {
				RHS[noE] *= (-1);
				flipFlg =1;
			    }
		    }
	    }//while(sti.hasMoreTokens())
	for(i=0;i<noCV;i++)
	    {
		for(j=0;j<noV;j++)
		    {
			if(coVar[i].endsWith(vars[j]))//Validate coeff*Var
			    {
				if(coVar[i].equals(vars[j]))
				    {
					if(flipFlg == 1)
					    sTable[noE][j] += -1;
					else
					    sTable[noE][j] += 1;
					incFlg = 1;
				    }
				else
				    {
					tmp = coVar[i].substring(0,(coVar[i].lastIndexOf(vars[j])));
					//System.out.println("tmp:"+tmp);
					if(tmp.equals("-"))
					    tmp+="1";
					if(flipFlg == 1)
					    sTable[noE][j] += (-1)*Double.parseDouble(tmp);
					else
					    sTable[noE][j] += Double.parseDouble(tmp);
					incFlg=1;
				    }
				break;
			    }
		    }
	    }
	if(incFlg == 1)
	    {
		if(flipFlg == 1)// mult by -1
		    {
			if(coVar[noCV].equals("="))
			    slk[noE] = 0;
			else if(coVar[noCV].equals("<="))
			    {slk[noE] = -1; TnoV++;}
			else if(coVar[noCV].equals(">="))
			    slk[noE] = 1; 
		    }
		else
		    {
			if(coVar[noCV].equals("="))
			    slk[noE] = 0;
			else if(coVar[noCV].equals("<="))
			    slk[noE] = 1;
			else if(coVar[noCV].equals(">="))
			    {slk[noE] = -1;TnoV++;}
		    }
		TnoV++;
		noE++;
	    }
    }
}
