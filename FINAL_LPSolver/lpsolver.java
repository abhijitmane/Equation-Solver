import java.io.*;
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.util.ArrayList;
import javax.swing.filechooser.*;
import java.awt.datatransfer.*;
import java.util.*;
import java.io.File;
import javax.swing.text.DefaultEditorKit;
import javax.swing.text.TextAction;

public class lpsolver extends JFrame implements ActionListener {
    addVar ad;
    simplex sm = new simplex();
    JMenuBar jmenubar = new JMenuBar();
    JFileChooser chooser = new JFileChooser();
    JTextPane p1 = new JTextPane();
    JTextPane p2 = new JTextPane();
    JTextPane p3 = new JTextPane();
    JTextPane p4 = new JTextPane();
    JTextPane p5 = new JTextPane();
    JMenu jmenu1 = new JMenu("File");
    JMenuItem jmenuitem1 = new JMenuItem("New");
    JMenuItem jmenuitem2 = new JMenuItem("Open");
    JMenuItem jmenuitem3 = new JMenuItem("Save");
    JMenuItem jmenuitem4 = new JMenuItem("Save As");
    JMenuItem jmenuitem5 = new JMenuItem("Exit");
    JMenu jmenu2 = new JMenu("Edit");
	
    JMenuItem jmenuitem6 = new JMenuItem("Cut");
    JMenuItem jmenuitem7 = new JMenuItem("Copy");
    JMenuItem jmenuitem8 = new JMenuItem("Paste");
    JMenuItem jmenuitem9 = new JMenuItem("Select All");
    JMenuItem jmenuitem10 = new JMenuItem("Find");
    JMenu jmenu3 = new JMenu("Solve");
    JMenuItem jmenuitem11 = new JMenuItem("Simplex");
    JMenuItem jmenuitem12 = new JMenuItem("Duality");
    JMenu jmenu6 = new JMenu("Sensitivity");
    JCheckBoxMenuItem jc1 = new JCheckBoxMenuItem("Range of Values for Objective Function Coefficients of NBVs");
    JCheckBoxMenuItem jc2 = new JCheckBoxMenuItem("Range of Values for Objective Function Coefficients of BVs");
    JCheckBoxMenuItem jc3 = new JCheckBoxMenuItem("Range of Values for RHS values of Constraints");
    JCheckBoxMenuItem jc4 = new JCheckBoxMenuItem("Range of Values for Constraint Coefficients of NBVs");
    JMenuItem jc5 = new JMenuItem("Add a New Decision Variable");
    JMenu jmenu4 = new JMenu("Help");
    JMenuItem jmenuitem14 = new JMenuItem("Tutorial");
    JMenuItem jmenuitem15 = new JMenuItem("Features");
    JMenu jmenu5 = new JMenu("About");
    JMenuItem jmenuitem16 = new JMenuItem("Author(s)");
    JMenuItem jmenuitem17 = new JMenuItem("Project");
    public Clipboard system;
    public StringSelection stsel;
    String rowstring,value;
    JTabbedPane jtp = new JTabbedPane(SwingConstants.BOTTOM);
    ArrayList list = new ArrayList();
    String path="new",s1="",s2="",s3="",s4="",s5="",sens="";
    int newFileFlg = 0;
    lpsolver() {
	setTitle("LPSolver-["+path+"]");
	Container contentpane = getContentPane();
	setSize(1300,800);
	addWindowListener(new WindowAdapter () { public void windowClosing (WindowEvent e) { System.exit(0); } } );
	jmenu1.add(jmenuitem1);
	jmenu1.add(jmenuitem2);
	jmenu1.add(jmenuitem3);
	jmenu1.add(jmenuitem4);
	jmenu1.addSeparator();
	jmenu1.add(jmenuitem5);	     
	jmenuitem1.addActionListener(this);
	jmenuitem2.addActionListener(this);
	jmenuitem3.addActionListener(this);
	jmenuitem4.addActionListener(this);
	jmenuitem5.addActionListener(this);
	jmenubar.add(jmenu1);

	jmenu2.add(jmenuitem6);
	jmenu2.add(jmenuitem7);
	jmenu2.add(jmenuitem8);
	jmenu2.add(jmenuitem9);
	jmenu2.add(jmenuitem10);
	jmenuitem6.addActionListener(this);
	jmenuitem7.addActionListener(this);
	jmenuitem8.addActionListener(this);
	jmenuitem9.addActionListener(this);
	jmenuitem10.addActionListener(this);
	jmenubar.add(jmenu2);
	
	jmenu3.add(jmenuitem11);
	jmenu3.add(jmenuitem12);
	jmenuitem11.addActionListener(this);
	jmenuitem12.addActionListener(this);
	jmenubar.add(jmenu3);
	
	jmenu6.add(jc1);
	jmenu6.add(jc2);
	jmenu6.add(jc3);
	jmenu6.add(jc4);
	jmenu6.add(jc5);
	jc1.addActionListener(this);
	jc2.addActionListener(this);
	jc3.addActionListener(this);
	jc4.addActionListener(this);
	jc5.addActionListener(this);
	jmenubar.add(jmenu6);

	jmenu4.add(jmenuitem14);
	jmenu4.add(jmenuitem15);
	jmenuitem14.addActionListener(this);
	jmenuitem15.addActionListener(this);
	jmenubar.add(jmenu4);

	jmenu5.add(jmenuitem16);
	jmenu5.add(jmenuitem17);
	jmenuitem16.addActionListener(this);
	jmenuitem17.addActionListener(this);
	jmenubar.add(jmenu5);
	
	jmenu1.setMnemonic(KeyEvent.VK_F);
	jmenu2.setMnemonic(KeyEvent.VK_E);
	jmenu3.setMnemonic(KeyEvent.VK_V);
	jmenu4.setMnemonic(KeyEvent.VK_H);
	jmenu5.setMnemonic(KeyEvent.VK_A);
	jmenu6.setMnemonic(KeyEvent.VK_S);
	jmenuitem1.setMnemonic(KeyEvent.VK_N);
	jmenuitem1.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N,ActionEvent.CTRL_MASK));
	jmenuitem2.setMnemonic(KeyEvent.VK_O);
	jmenuitem2.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O,ActionEvent.CTRL_MASK));
	jmenuitem3.setMnemonic(KeyEvent.VK_S);
	jmenuitem3.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,ActionEvent.CTRL_MASK));
	jmenuitem4.setMnemonic(KeyEvent.VK_A);
	jmenuitem4.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A,ActionEvent.ALT_MASK));
	jmenuitem5.setMnemonic(KeyEvent.VK_X);
	jmenuitem5.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X,ActionEvent.ALT_MASK));
	jmenuitem6.setMnemonic(KeyEvent.VK_X);
	jmenuitem6.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X,ActionEvent.CTRL_MASK));
	jmenuitem7.setMnemonic(KeyEvent.VK_C);
	jmenuitem7.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C,ActionEvent.CTRL_MASK));
	jmenuitem8.setMnemonic(KeyEvent.VK_V);
	jmenuitem8.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V,ActionEvent.CTRL_MASK));
	jmenuitem9.setMnemonic(KeyEvent.VK_A);
	jmenuitem9.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A,ActionEvent.CTRL_MASK));
	jmenuitem10.setMnemonic(KeyEvent.VK_F);
	jmenuitem10.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F,ActionEvent.CTRL_MASK));
	jmenuitem11.setMnemonic(KeyEvent.VK_X);
	jmenuitem11.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X,ActionEvent.ALT_MASK));
	jmenuitem12.setMnemonic(KeyEvent.VK_D);
	jmenuitem12.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D,ActionEvent.ALT_MASK));
	jc5.setMnemonic(KeyEvent.VK_A);
	jc5.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A,ActionEvent.ALT_MASK));
	jmenuitem14.setMnemonic(KeyEvent.VK_T);
	jmenuitem14.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T,ActionEvent.ALT_MASK));
	jmenuitem15.setMnemonic(KeyEvent.VK_F);
	jmenuitem15.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F,ActionEvent.ALT_MASK));
	jmenuitem16.setMnemonic(KeyEvent.VK_A);
	jmenuitem16.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A,ActionEvent.ALT_MASK));
	jmenuitem17.setMnemonic(KeyEvent.VK_P);
	jmenuitem17.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_P,ActionEvent.ALT_MASK));
			
	setJMenuBar(jmenubar);
	jtp.addTab("Problem", new ProbPanel(p1));
	jtp.addTab("Simplex", new SimPanel(p2));
	jtp.addTab("Dual", new DualPanel(p3));
	jtp.addTab("Sensitivity", new SensPanel(p4));
	jtp.addTab("Errors", new ErrPanel(p5));
	contentpane.add(jtp);
	setVisible(true);
    }
    public void actionPerformed(ActionEvent e) {
	int result;
	JMenuItem jmenuitem = (JMenuItem)e.getSource();
	String str = jmenuitem.getText();
	JFileChooser chooser = new JFileChooser();
	if(str == "New"){
	    p1.setText("");
	    p2.setText("");
	    p3.setText("");
	    p4.setText("");
	    p5.setText("");
	    sm.simpSolveFlg = 0;
	    newFileFlg = 0;
	    s1="";s2="";s3="";s4="";s5="";sens="";
	    path = "new";
	    setTitle("LPSolver-["+path+"]");
	}
	if(str == "Open"){
	    result=chooser.showOpenDialog(null);
	    if(result==JFileChooser.APPROVE_OPTION)
                {
		    try{
			String data ="";
			File f1=chooser.getSelectedFile();
			path=f1.getAbsolutePath();
			FileInputStream fi = new FileInputStream(path);
			f1 = new File(path);
			long i = f1.length();
			byte b[] = new byte[(int)i];
			fi.read(b);
			data+= new String(b,0,(int)i);
			fi.close();
			p1.setText(data);
			p2.setText("");
			p3.setText("");
			p4.setText("");
			p5.setText("");
			sm.simpSolveFlg = 0;
			newFileFlg = 1;
			s1="";s2="";s3="";s4="";s5="";sens="";
			setTitle("LPSolver-["+path+"]");
		    }catch(Exception ex){}
		}
	}
    	if(str == "Save"){
	    try{
		String data="";
		if(newFileFlg == 0)
		    {
			newFileFlg = 1;
			result=chooser.showSaveDialog(null);
			if(result==JFileChooser.APPROVE_OPTION)
			    {
				File f1=chooser.getSelectedFile();
				path=f1.getAbsolutePath();
			    }
		    }
		if(newFileFlg == 1)
		    {
			File f = new File(path);
			if(f.exists())
			    f.delete();
			RandomAccessFile fp =new RandomAccessFile(path,"rw");
			data += p1.getText();
			fp.writeBytes(data);
			fp.close();
			setTitle("LPSolver-["+path+"]");
		    }
	    }
	    catch(Exception ex){}
	}
	if(str == "Save As"){
	    result=chooser.showSaveDialog(null);
            if(result==JFileChooser.APPROVE_OPTION)
		{
		    try{
			File f1=chooser.getSelectedFile();
			path=f1.getAbsolutePath();
			RandomAccessFile fp =new RandomAccessFile(path,"rw");
			String data = p1.getText();
			fp.writeBytes(data);
			fp.close();
			setTitle("LPSolver-["+path+"]");
		    }
		    catch(Exception ex){}
		}
	}
	if(str == "Exit" ){
	    System.exit(0);
	}
	
	if(str == "Cut"){
	    cut();  
	}
	if(str == "Copy"){
	    copy();
	}
	if(str == "Paste"){
	    paste();
	}
	if(str == "Select All"){
	    selectAll();
	}
	if(str == "Find"){
	    search();
	}
	if(str == "Simplex"){
	    String input = p1.getText();
	    String fname = "simplex.txt";
	    try
		{
		    File f1 = new File(fname);
		    if(f1.exists())
			f1.delete();
		    RandomAccessFile fp =new RandomAccessFile(fname,"rw");
		    fp.writeBytes(input);
		    fp.close();
		    s1="";s2="";s3="";s4="";s5="";sens="";
		    input = sm.callSimplex(fname,"simplex",sm);
		    if((input.substring(0,7)).equals("ERRORS:")){
			p2.setText("");
			p5.setText(input);
		    }
		    else
			{
			    p2.setText(input);
			    p5.setText("");
			    sm.simpSolveFlg = 1;
			}
		}
	    catch(Exception ex){}
	}
	if(str == "Duality"){
	    String input = p1.getText();
	    String fname = "simplex.txt";
	    try
		{
		    File f1 = new File(fname);
		    if(f1.exists())
			f1.delete();
		    RandomAccessFile fp =new RandomAccessFile(fname,"rw");
		    fp.writeBytes(input);
		    fp.close();
		    
		    input = sm.callSimplex(fname,"dual",sm);
		    if((input.substring(0,7)).equals("ERRORS:")){
			p3.setText("");
			p5.setText(input);
		    }
		    else
			{
			    p3.setText(input);
			    p5.setText("");
			}
		}
	    catch(Exception ex){}
	}
	
	if(str == "Tutorial"){
	    try{
		java.io.File  opChm = new java.io.File("lp.chm");
		Desktop.getDesktop().open(opChm.getAbsoluteFile());
	    }
	    catch(Exception ex){}
	    

	    
	    /*
		FileInputStream fi = new FileInputStream("help.txt");
		File f1 = new File("help.txt");
		long i = f1.length();
		byte b[] = new byte[(int)i];
		fi.read(b);
		data+= new String(b,0,(int)i);
		fi.close();
		JOptionPane.showMessageDialog(null,"<html><body color=black><font size=4><b>"+data+"","Help Tutorial",JOptionPane.INFORMATION_MESSAGE);
	    }
	    catch(Exception ex){}*/
	}
	if(str == "Features"){
	    try{
		java.io.File  opChm = new java.io.File("lp.chm");
		Desktop.getDesktop().open(opChm.getAbsoluteFile());
	    }
	    catch(Exception ex){}

	    
	}
	if(str == "Author(s)"){
	    JOptionPane.showMessageDialog(null,"<html><body color=black><font size=4><b>08231:Sagar Takawane\n<html><body color=black><font size=4><b>07123:Abhijit Mane\n<html><body color=black><font size=4><b>08123:Rajesh Kumar","About Us",JOptionPane.INFORMATION_MESSAGE);
	}
	if(str == "Project"){
	    try{
		java.io.File  opChm = new java.io.File("help.html");
		Desktop.getDesktop().open(opChm.getAbsoluteFile());
	    }
	    catch(Exception ex){}
	}
	if(jmenuitem == jc1){ // in Obj, coeff. of NBV 
	    if(jc1.getState() == true){
		if((sm.iFlg == 1) || (sm.uFlg == 1) || (sm.cFlg == 1))
		    {
			if(sm.iFlg == 1)
			    {
				s1 = "Error!!! see Error Tab";
				p5.setText("ERROR:\n\tSolution is Infeasible, So unable to do Sensitivity analysis\n ");
			    }
			if(sm.uFlg == 1)
			    {
				s1 = "Error!!! see Error Tab";
				p5.setText("ERROR:\n\tSolution is Unbounded, So unable to do Sensitivity analysis\n");
			    }
			if(sm.cFlg == 1)
			    {
				s1 = "Error!!! see Error Tab";
				p5.setText("ERROR:\n\tCycle occurs... No solution, So unable to do Sensitivity analysis\n");
			    }
		    }
		else if(sm.simpSolveFlg == 0)
		    p5.setText("ERROR:\n\tProblem is not solved using Simplex Method\n");
		else
		    s1 = sm.objNBV();
	    }
	    else{
		s1 = "";
	    }
	    sens = s1 + s2 + s3 + s4 + s5;
	    p4.setText(sens);
	    }
	if(jmenuitem == jc2){ // in Obj, coeff. of BV 
	    if(jc2.getState() == true){
		if((sm.iFlg == 1) || (sm.uFlg == 1) || (sm.cFlg == 1))
		    {
			if(sm.iFlg == 1)
			    {
				s2 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Infeasible, So unable to do Sensitivity analysis\n ");
			    }
			if(sm.uFlg == 1)
			    {
				s2 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Unbounded, So unable to do Sensitivity analysis\n");
			    }
			if(sm.cFlg == 1)
			    {
				s2 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tCycle occurs... No solution, So unable to do Sensitivity analysis\n");
			    }
		    }
		else if(sm.simpSolveFlg == 0)
		    {
			s2 = "Error!!! see Error Tab\n";
			p5.setText("ERROR:\n\tProblem is not solved using Simplex Method\n");
			JOptionPane.showMessageDialog(null,"<html><body color=black><font size=4><b>First Solve the Problem,\n then go for Sensitivity Analysis","Warning!!!",JOptionPane.INFORMATION_MESSAGE);
		    }
		else
		    s2 = sm.objBV();
	    }
	    else{
		s2 = "";
	    }
	    sens = s1 + s2 + s3 + s4 + s5;
	    p4.setText(sens);
	}
	if(jmenuitem == jc3){ // RHS of constraint
	    if(jc3.getState() == true){
		if((sm.iFlg == 1) || (sm.uFlg == 1) || (sm.cFlg == 1))
		    {
			if(sm.iFlg == 1)
			    {
				s3 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Infeasible, So unable to do Sensitivity analysis\n ");
			    }
			if(sm.uFlg == 1)
			    {
				s3 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Unbounded, So unable to do Sensitivity analysis\n");
			    }
			if(sm.cFlg == 1)
			    {
				s3 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tCycle occurs... No solution, So unable to do Sensitivity analysis\n");
			    }
		    }
		else if(sm.simpSolveFlg == 0)
		    {
			s3 = "Error!!! see Error Tab\n";
			p5.setText("ERROR:\n\tProblem is not solved using Simplex Method\n");
		    }
		else
		    s3 = sm.cnstrRHS();
	    }
		else{
		    s3 = "";
		}
		sens = s1 + s2 + s3 + s4 + s5;
		p4.setText(sens);    
	    }
	    if(jmenuitem == jc4){ //in constraint, coeff. of NBV 
		if(jc4.getState() == true){
		    if((sm.iFlg == 1) || (sm.uFlg == 1) || (sm.cFlg == 1))
		    {
			if(sm.iFlg == 1)
			    {
				s4 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Infeasible, So unable to do Sensitivity analysis\n ");
			    }
			if(sm.uFlg == 1)
			    {
				s4 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Unbounded, So unable to do Sensitivity analysis\n");
			    }
			if(sm.cFlg == 1)
			    {
				s4 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tCycle occurs... No solution, So unable to do Sensitivity analysis\n");
			    }
		    }
		    else if(sm.simpSolveFlg == 0)
			{
			    s4 = "Error!!! see Error Tab\n";
			    p5.setText("ERROR:\n\tProblem is not solved using Simplex Method\n");
			}
		    else
			s4 = sm.cnstrNBV();
		}
		else{
		    s4 = "";
		}
		sens = s1 + s2 + s3 + s4 + s5;
		p4.setText(sens);
	    }
	    if(jmenuitem == jc5){ // add new var
		if((sm.iFlg == 1) || (sm.uFlg == 1) || (sm.cFlg == 1))
		    {
			if(sm.iFlg == 1)
			    {
				s5 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Infeasible, So unable to do Sensitivity analysis\n ");
			    }
			if(sm.uFlg == 1)
			    {
				s5 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tSolution is Unbounded, So unable to do Sensitivity analysis\n");
			    }
			if(sm.cFlg == 1)
			    {
				s5 = "Error!!! see Error Tab\n";
				p5.setText("ERROR:\n\tCycle occurs... No solution, So unable to do Sensitivity analysis\n");
			    }
		    }
		else if(sm.simpSolveFlg == 0)
		    {
			p5.setText("ERROR:\n\tProblem is not solved using Simplex Method\n");
			s5 = "Error!!! see Error Tab\n";
		    }
		else
		    ad = new addVar(sm,this);
		
	    }
    }
    
    public static void main(String args[]) {
	new lpsolver();
    }

    public void cut()
    {
        p1.cut();
    }

    public void copy()
    {
        p1.copy();
    }

    public void paste()
    {
        p1.paste();
    }
    public void selectAll()
    {
	p1.selectAll();
    }
    public void search()
    {
	int index,end;
	JFrame frame = new JFrame();
        p1.requestFocus();
        String item = JOptionPane.showInputDialog("Enter Search String");
        String text = p1.getText();
        index = text.indexOf(item);
        end = index + item.length();
        p1.select(index, end);
        if(index == -1)
            JOptionPane.showMessageDialog(frame, "String not found!");
    
    }
    
    public void save() 
    {
	int result=chooser.showSaveDialog(null);
	if(result==JFileChooser.APPROVE_OPTION)
	    {
		try{
		    String data;
		    File f1=chooser.getSelectedFile();
		    String path=f1.getAbsolutePath();
		    File f = new File(path);
		    if(f.exists())
			f.delete();
		    
		    RandomAccessFile fp =new RandomAccessFile(path,"rw");
		    data = "PROBLEM:\n";
		    data += p1.getText();
		    data += "\n\nSIMPLEX:\n";
		    data += p2.getText();
		    data += "\n\nDUAL:\n";
		    data += p3.getText();
		    data += "\n\nSENSITIVITY:\n";
		    data += p4.getText();
		    data += "\n\nERRORS:\n";
		    data += p5.getText();
		    fp.writeBytes(data);
		    fp.close();
		}
		catch(Exception ex){}
	    }

    }
    
}

class ProbPanel extends JPanel 
{
    public ProbPanel(JTextPane p1) 
    {
	setLayout(new BorderLayout());
	
	JScrollPane jsp=new JScrollPane(p1,ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
	add(jsp);
	p1.setFont(new Font("Serif",Font.BOLD,15));
    }
}


class SimPanel extends JPanel 
{
    public SimPanel(JTextPane p2) 
    {
	setLayout(new BorderLayout());
	JScrollPane jsp=new JScrollPane(p2,ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
	add(jsp);
	p2.setFont(new Font("Serif",Font.BOLD,15));
	p2.setEditable(false);
    }
}

class DualPanel extends JPanel 
{
    public DualPanel(JTextPane p3) 
    {
	setLayout(new BorderLayout());
	JScrollPane jsp=new JScrollPane(p3,ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
	add(jsp);
	p3.setFont(new Font("Serif",Font.BOLD,15));
	p3.setEditable(false);
    }
}

class SensPanel extends JPanel 
{
    public SensPanel(JTextPane p4) 
    {
	setLayout(new BorderLayout());
	JScrollPane jsp=new JScrollPane(p4,ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
	add(jsp);
	p4.setFont(new Font("Serif",Font.BOLD,15));
	p4.setEditable(false);
    }
}

class ErrPanel extends JPanel 
{
    public ErrPanel(JTextPane p5) 
    {
	setLayout(new BorderLayout());
	JScrollPane jsp=new JScrollPane(p5,ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
	add(jsp);
	p5.setFont(new Font("Serif",Font.BOLD,15));
	p5.setEditable(false);
    }
}

