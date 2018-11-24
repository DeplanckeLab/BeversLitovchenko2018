package mitorepeats;

/**
 *
 * @author litovche
 */
public class Repeat {
    Integer ATs = 0;
    String Spacer = "";
    Integer TAs = 0;
    
    String repeatCollapsed = "";
    
    public Repeat () {}
    
    public Repeat (String s) {
        int i = 3; // because CTA is always first 3 letters
        while (s.substring(i, i + 2).equals("AT") &
               !s.substring(i + 1, i + 2).equals("G")) {
            ATs++;
            i = i + 2;
        }
        while (!s.substring(i, i + 2).equals("TA") &
               !s.substring(i + 1, i + 2).equals("G")) {
            Spacer = Spacer + s.substring(i, i + 1);
            i = i + 1;
        }
        while (s.substring(i, i + 2).equals("TA") &
               !s.substring(i + 1, i + 2).equals("G")) {
            TAs++;
            i = i + 2;
        }
        CollapseRepeat();
    }
    
    private void CollapseRepeat() {
        this.repeatCollapsed = "(AT)" + ATs + "_" + Spacer + "_(TA)" + TAs;
    }
}
