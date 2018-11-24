package mitorepeats;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author litovche
 */
public class Main {
    /**
     * @param filePath
     * @return 
     */
    
    private static HashMap sortByValues(HashMap map) { 
       List list = new LinkedList(map.entrySet());
       // Defined Custom Comparator here
       Collections.sort(list, new Comparator() {
            public int compare(Object o1, Object o2) {
                int revertResult = ((Comparable) ((Map.Entry) (o1)).getValue())
                  .compareTo(((Map.Entry) (o2)).getValue());
               if (revertResult == -1) {
                   return 1;
               }  else {
                   if (revertResult == 1) {
                       return -1;
                   } else {
                       return 0;
                   }
               }
            }
       });

       // Here I am copying the sorted list in HashMap
       // using LinkedHashMap to preserve the insertion order
       HashMap sortedHashMap = new LinkedHashMap();
       for (Iterator it = list.iterator(); it.hasNext();) {
              Map.Entry entry = (Map.Entry) it.next();
              sortedHashMap.put(entry.getKey(), entry.getValue());
       } 
       return sortedHashMap;
  }

    
    public static Map readFile (File filePath) {
        String line = null;
        Map<String, Integer> allRepeats = new HashMap<>();
        
        try { // FileReader reads text files in the default encoding.
            FileReader fileReader = new FileReader(filePath);

            // Always wrap FileReader in BufferedReader.
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            //System.out.println(filePath.getName());
            while((line = bufferedReader.readLine()) != null) {
                String[] data = line.split(" ");
                Repeat rep = new Repeat(data[1]);
                if (allRepeats.containsKey(rep.repeatCollapsed)) {
                    allRepeats.put(rep.repeatCollapsed, 
                                   allRepeats.get(rep.repeatCollapsed) + 
                                   Integer.parseInt(data[0]));
                } else {
                    allRepeats.put(rep.repeatCollapsed, Integer.parseInt(data[0]));
                }
                
                /*System.out.println(line);
                System.out.println(Integer.parseInt(data[0]) + " " + data[1] + 
                                   " " + rep.repeatCollapsed);*/
            }
            // Always close files.
            bufferedReader.close();         
        }
        catch(FileNotFoundException ex) {
            System.out.println("Unable to open file '" + filePath + "'");                
        }
        catch(IOException ex) {
            System.out.println("Error reading file '" + filePath + "'");    
        }
        return allRepeats;
    }
    
    public static void main(String[] args) {
        File folder = new File("/home/litovche/Desktop/testRepeats");
        File[] listOfFiles = folder.listFiles();
        
        ArrayList<String> allTheRepeats = new ArrayList<>();
        Map<String, Map<String, Integer>> fileStats = new HashMap<>();
        
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
                Map<String, Integer> repeatsInFile = readFile(listOfFiles[i]);
                /*System.out.println(listOfFiles[i].getName());
                for (String s : repeatsInFile.keySet()) {
                    System.out.println(s + " " + repeatsInFile.get(s));
                }*/
                allTheRepeats.addAll(repeatsInFile.keySet());
                fileStats.put(listOfFiles[i].getName(), repeatsInFile);
            }
        }
        
        Set<String> uniqueRepeats = new HashSet<String>(allTheRepeats);
        HashMap<String, Integer> occurrences = new HashMap();
        
        for (Iterator<String> it = uniqueRepeats.iterator(); it.hasNext(); ) {
            String f = it.next();
            occurrences.put(f, Collections.frequency(allTheRepeats, f));
        }
                      
        occurrences = sortByValues(occurrences);
        
        System.out.print("File\t");
        for (String key : occurrences.keySet()) {
            System.out.print(key + "\t");
        }
        System.out.println();
        for (String fileName : fileStats.keySet()) {
            System.out.print(fileName + "\t");
            Map<String, Integer> oneFileStats = fileStats.get(fileName);
            for (String key : occurrences.keySet()) {
                if (oneFileStats.containsKey(key)) {
                    System.out.print(oneFileStats.get(key) + "\t");
                } else {
                    System.out.print("0\t");
                }
            }
            System.out.println("\t");
        }
    }    
}
