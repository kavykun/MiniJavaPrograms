import java.util.*;
public class Palindrome {
  /**
   *
   */ 
  public static void main(String [] args) {
    userInterface();

  }
  /**
   *
   */
  public static void userInterface(){
    Scanner in = new Scanner(System.in);
    System.out.print("Palindrome?: ");
    String s = in.nextLine();
    String ls = s.toLowerCase();
    String sr = reverse(s);
    String lsr = sr.toLowerCase();
    int count = 0;
    for (int i = 0; i < ls.length(); i++){
      char c1 = ls.charAt(i);
      char c2 = lsr.charAt(i);
      if (c1 == c2) {
        count++ ;
      } 
    }
    if (count == ls.length()){
      System.out.println(s + " is a Palindrome.");
    } else {
        System.out.println(s + " is not a Palindrome");
    }
    System.out.println(ls);
    System.out.println(lsr);
  }
  /**
   *
   */ 
  public static String reverse(String s){
    String p = "";
    for (int i = s.length(); i > 0; i--){
      p += s.charAt(i-1);
    }
    return p;

}

}

