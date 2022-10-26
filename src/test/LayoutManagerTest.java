package test;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import javax.swing.JButton;
import javax.swing.JFrame;

public class LayoutManagerTest {

    public static void main( String[] args ) {

        JFrame f1 = new JFrame( "BorderLayout" );
        f1.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );

        f1.add( new JButton( "btn1" ) );
        f1.add( new JButton( "btn2" ) );
        f1.add( new JButton( "btn3" ) );
        f1.add( new JButton( "btn4" ) );
        f1.add( new JButton( "btn5" ) );
        f1.setSize( 500, 200 );
        f1.setLocationRelativeTo( null );



        JFrame f2 = new JFrame( "BorderLayout with regions" );
        f2.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );

        f2.add( new JButton( "btn1" ), BorderLayout.NORTH );
        f2.add( new JButton( "btn2" ), BorderLayout.SOUTH );
        f2.add( new JButton( "btn3" ), BorderLayout.WEST );
        f2.add( new JButton( "btn4" ), BorderLayout.EAST );
        f2.add( new JButton( "btn5" ), BorderLayout.CENTER );
        f2.setSize( 500, 200 );
        f2.setLocationRelativeTo( null );



        JFrame f3 = new JFrame( "FlowLayout" );
        f3.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );
        f3.setLayout( new FlowLayout() );

        f3.add( new JButton( "btn1" ) );
        f3.add( new JButton( "btn2" ) );
        f3.add( new JButton( "btn3" ) );
        f3.add( new JButton( "btn4" ) );
        f3.add( new JButton( "btn5" ) );
        f3.setSize( 500, 200 );
        f3.setLocationRelativeTo( null );



        f1.setVisible( true );
        f2.setVisible( true );
        f3.setVisible( true );

    }

}