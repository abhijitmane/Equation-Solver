import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
public class addVar extends JFrame implements ActionListener 
{
    simplex sm;
    lpsolver lp;
    String sRet,vname,objf,cnstr,restr="0";
    public static addVar instance;
    public static JLabel label1;
    public static JLabel label2;
    public static JLabel label3;
    public static JLabel label4;
    public static JTextField text5;
    public static JTextField text6;
    public static JTextField text7;
    public static JRadioButton box8;
    public static JRadioButton box9;
    public static JButton button10;
    public static JButton button11;
    public ButtonGroup bg = new ButtonGroup();
    public addVar(simplex sm,lpsolver lp)
    {
	this.sm = sm;
	this.lp = lp;
	setTitle("Add a New Variable");
	setSize(550,200);
	setLayout(null);
	addWindowListener(new WindowAdapter () { public void windowClosing (WindowEvent e) { setVisible(false); } } );
	label1 = new JLabel("Variable Name:");
	add(label1);
	label1.setBounds(7,22,120,20);
	label2 = new JLabel("Objective Function Value:");
	add(label2);
	label2.setBounds(6,47,200,20);
	label3 = new JLabel("Constraint Values:");
	add(label3);
	label3.setBounds(7,74,150,20);
	label4 = new JLabel("Constraint Values are entered in the order of constraints as a comma separated list");
	add(label4);
	label4.setBounds(7,101,530,20);
	text5 = new JTextField();
	add(text5);
	text5.setBounds(194,22,70,20);
	text6 = new JTextField();
	add(text6);
	text6.setBounds(193,49,70,20);
	text7 = new JTextField();
	add(text7);
	text7.setBounds(193,76,120,20);
	box8 = new JRadioButton(" >= 0");
	add(box8);
	box8.setBounds(17,134,70,20);
	box8.addActionListener(this);
	box9 = new JRadioButton(" <= 0");
	add(box9);
	box9.setBounds(111,134,70,20);
	box9.addActionListener(this);
	button10 = new JButton("Ok");
	button10.addActionListener(this);
	add(button10);
	button10.setBounds(257,133,80,23);
	button11 = new JButton("Cancel");
	button11.addActionListener(this);
	add(button11);
	button11.setBounds(377,133,80,23);
	bg.add(box8);
	bg.add(box9);
	setVisible(true);
    }
    
    public void actionPerformed (ActionEvent e) {
	if(e.getSource() == button10) {
	    sRet = "";
	    vname = text5.getText();
	    objf = text6.getText();
	    cnstr = text7.getText();
	    sRet = sm.addNewVar(this,sm);
	    if(sRet.indexOf("No Change in Solution")!= -1)
		lp.p4.setText(lp.s1+lp.s2+lp.s3+lp.s4+"\n\n* Add New Variable:-\n"+sRet); 
	    else
		{
		    lp.p4.setText(lp.s1+lp.s2+lp.s3+lp.s4+"\n\n* Variable "+vname+" is added successfully, see updated solution on SIMPLEX tab\n"); 
		    lp.p2.setText(lp.p2.getText()+"\n\n* Add New Variable:-\n"+sRet); 
		}
	    lp.p5.setText("");
	    setVisible(false);
	}
	if(e.getSource() == button11) {
	    setVisible(false);
	}
	if(e.getSource() == box8) {
	    restr = "1";
	}
	if(e.getSource() == box9) {
	    restr = "-1";
	}
    }
    
    /*public static void main (String [] args) 
    {
	instance = new addVar();
	//instance.show();
	}*/
}
