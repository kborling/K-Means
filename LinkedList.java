    /*
   * Kevin Borling
   * CSCD 429
   * LinkedList and Node class to load data and test set into.
   *
   */
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
           
   public class LinkedList
      {
      private Node head;
      private int size;
      
   //--------------------------Node--------------------------//
      public class Node
         {

         protected String data;
         protected String geneID;
         protected String essential;
         protected String classify;
         protected String complex;
         protected String phenotype;
         protected String motif;
         protected String chromosome;
         protected String function;
         protected String localization;
         protected float weight;
         protected Node next;
         
         protected Node(String geneID, String essential, String classify, String complex, String phenotype, String motif, String chromosome, String function, String localization, Node next)
         {
            this.geneID = geneID;
            this.essential = essential;
            this.classify = classify;
            this.complex = complex;
            this.phenotype = phenotype;
            this.motif = motif;
            this.chromosome = chromosome;
            this.function = function;
            this.localization = localization;
            this.weight = 0.0f;
            this.next = next;
         }// End EVC
            
         protected Node(String geneID, String essential, String classify, String complex, String phenotype, String motif, String chromosome, String function, String localization)
         {
            this(geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization, null);
         }// End DVC
      
         public void setNext(Node next) 
         { 
            this.next = next;
         }// End setNext
         public Node getNext()
         {
            return this.next;
         }// End getNext

         public String getGeneID()
         {
            return this.geneID;
         }// End getGeneID
         public String getEssential()
         {
            return this.essential;
         }// End getEssential
         public String getClassify()
         {
            return this.classify;
         }// End getClassify
         public String getComplex()
         {
            return this.complex;
         }// End getComplex
         public String getPhenotype()
         {
            return this.phenotype;
         }// End getPhenotype
         public String getMotif()
         {
            return this.motif;
         }// End getMotif
         public String getChromosome()
         {
            return this.chromosome;
         }// End getChromosome
         public String getFunction()
         {
            return this.function;
         }// End getFunction
         public String getLocalization()
         {
            return this.localization;
         }// End getLocalization
   
         public void setGeneID(String geneID)
         {
            this.geneID = geneID;
         }// End setGeneID
         public void setEssential(String essential)
         {
            this.essential = essential;
         }// End setEssential
         public void setClassify(String classify)
         {
            this.classify = classify;
         }// End setClassify
         public void setComplex(String complex)
         {
            this.complex = complex;
         }// End setComplex
         public void setPhenotype(String phenotype)
         {
            this.phenotype = phenotype;
         }// End setPhenotype
         public void setMotif(String motif)
         {
            this.motif = motif;
         }// End setMotif
         public void setChromosome(String chromosome)
         {
            this.chromosome = chromosome;
         }// End setChromosome
         public void setFunction(String function)
         {
            this.function = function;
         }// End setFunction
         public void setLocalization(String localization)
         {
            this.localization = localization;
         }// End setLocalization
      
      }// End Class




         //--------------------------LinkedList-----------------------//
      public LinkedList()
      {
         this.head = null;
         this.size = 0;
      }//End DVC
   		
      public LinkedList(int size)
      {
         this.size = size;
      }// End EVC

      public Node getHead()
      {
         return this.head;
      }
         
      public boolean headless()
      {
         return this.head == null;
      }// End isEmpty

      public void addHead(String geneID, String essential, String classify, String complex, String phenotype, String motif, String chromosome, String function, String localization)
      {
         Node myNode = new Node(geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization);
         myNode.next = this.head;
         this.head = myNode;
         this.size++;
      }//end addFirst method
   
      public void add (String geneID, String essential, String classify, String complex, String phenotype, String motif, String chromosome, String function, String localization)
      {
         if (headless())// Checks if head is null or empty
            addHead(geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization);
         else
            {
            Node temp = new Node(geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization);
            Node current = head;
         // Start at head node, move to the end of the list
            while(current.getNext() != null)
               {
               current = current.getNext();
               }
         // Sets last node's "next" reference to the new node
            current.setNext(temp);
            size++;
            }
         }// End add

      public void sort()
      {
         Node first, smallest, curr;
         String temp;
         
         for (first = this.head; first.next != null; first = first.next)
            {
               smallest = first; 
               for (curr = first.next; curr != null; curr = curr.next) 
                  if (curr.geneID.compareTo(smallest.geneID) < 0)
                     smallest = curr;
                  
               temp = first.geneID;
               first.geneID = smallest.geneID;
               smallest.geneID = temp;
            }
      }//End sort

      public void getEuclid(LinkedList data, LinkedList test)
      {
         Node dataHead = data.getHead();
         Node testHead = test.getHead();
         String local = "", gene = "";
         String one = "", two = "";
         float lowestDistance = 999, distance = 0.0f;
         int innercount = 0, outercount = 0;

         for (Node curr = testHead; curr != null; curr = curr.next)
         {
            outercount ++;
            one = oneLine(curr);

            for (Node search = dataHead; search != null; search = search.next)
            {
               innercount++;
               two = oneLine2(search);

              distance = calcEuclid(one, two);

              // Keep track of shortest Distance
              if(distance < lowestDistance)
              {
                  lowestDistance = distance;
                  local = search.localization;
                  gene = search.geneID;
              }
            }// End data set inner for

            // Sets localization of best match.
            curr.localization = local;
            curr.weight = lowestDistance;

            // Reset Values
            distance = 0.0f;
            lowestDistance = 999;
            local = "";
            gene = "";
         }// End test set outer for
      }// End calcEuclid


    public float calcEuclid(final String testString, final String dataString) 
    {


      // Parse File into tokens, splitting by commas or double quotations
      String testData, dataData;

       testData = testString.replace(".","");
       dataData = dataString.replace(".","");
      // Parse File into tokens, splitting by commas or double quotations
       String [] tokens1 = testData.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)");
       String [] tokens2 = dataData.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)");

        final ArrayList<String> testTokens = new ArrayList<String>();
        final ArrayList<String> dataTokens = new ArrayList<String>();

             for(String token1 : tokens1)
             {
               testTokens.add(token1);

             }
             for(String token2 : tokens2)
             {
               dataTokens.add(token2);
             }

        final Set<String> tokenSet = new HashSet<String>();

        tokenSet.addAll(testTokens);
        tokenSet.addAll(dataTokens);

        float distance = 0.0f;
        for (final String token : tokenSet) 
        {
            int testCount = 0;
            int dataCount = 0;
            for (final String sToken : testTokens) 
            {
                if (sToken.equals(token) && !sToken.equals("?")) // Ignore Mising Values
                    testCount++;
            }
            for (final String sToken : dataTokens) 
            {
                if (sToken.equals(token) && !sToken.equals("?"))
                    dataCount++;
            }

            distance += ((testCount - dataCount) * (testCount - dataCount));
        }

        return (float) Math.sqrt(distance);
    }// End 

      public String oneLine(Node curr)
      {
         String result = "";
            result += curr.geneID + "," + curr.essential + "," + curr.classify + "," + curr.complex + "," + curr.phenotype + "," + curr.motif + "," + curr.chromosome + "\n";
            return result;
      }// End oneLine
      public String oneLine2(Node curr)
      {
         String result = "";
            result += curr.geneID + "," + curr.essential + "," + curr.classify + "," + curr.complex + "," + curr.phenotype + "," + curr.motif + "," + curr.chromosome + "," + curr.function + "," + curr.localization + "\n";
            return result;
      }// End oneLine2

       public void removeDuplicates() 
       {

           Node prev = this.head;    
           Node curr = this.head.next;

           while(curr != null)
           {
               if(curr.geneID.equals(prev.geneID))
               {
                   prev.next = curr.next;
                   curr = curr.next;

               }else
               {
                   prev = curr;
                   curr = curr.next; 
               }
           }// End while
       }// End removeDuplicates


      @Override
      public String toString()
      {
      String result = "\n";
      for (Node curr = this.head; curr != null; curr = curr.next)
         result += curr.geneID + "," + curr.localization + "\n";
      return result;
      }//End toString
   }// End Class