using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using PythonPlotter;
using static System.Linq.Enumerable;
using CustomExtensions;

namespace numpy_functns
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            double value1 = 123.456;
            double value2 = 123.987;
            //
            // Take floors of these values.
            //
            double floor1 = Math.Floor(value1);
            double floor2 = Math.Floor(value2);
            double ceiling1 = Math.Ceiling(value1);

            //
            // Write first value and floor.
            //
            Console.WriteLine(value1);
            Console.WriteLine(floor1);
            Console.WriteLine(ceiling1);
            Console.WriteLine(value2);
            Console.WriteLine(floor2);
                               
            // Get ceiling of decimal value.
            decimal value3 = 456.789M;
            decimal ceiling2 = Math.Ceiling(value3);

            // Get ceiling of negative value.
            double value4 = -100.5;
            double ceiling3 = Math.Ceiling(value4);

            Console.WriteLine(value3);
            Console.WriteLine(ceiling2);
            Console.WriteLine(value4);
            Console.WriteLine(ceiling3);
        }

        private void button2_Click(object sender, EventArgs e)
        {
            // ... Create 2D array of strings.
            string[,] array = new string[,]
            {
            {"cat", "dog"},
            {"bird", "fish"},
            };
            // ... Print out values.
            Console.WriteLine(array[0, 0]);
            Console.WriteLine(array[0, 1]);
            Console.WriteLine(array[1, 0]);
            Console.WriteLine(array[1, 1]);

            // ... A one-dimensional array.
            int[] one = new int[2];
            one[0] = 1;
            one[1] = 2;
            Handle(one);

            // ... A two-dimensional array.
            int[,] two = new int[2, 2];
            two[0, 0] = 0;
            two[1, 0] = 1;
            two[0, 1] = 2;
            two[1, 1] = 3;
            Handle(two);
        }

    new static void Handle(Array array)
    {
        Console.WriteLine("Rank: " + array.Rank);
        switch (array.Rank)
        {
            case 1:
                for (int i = 0; i < array.Length; i++)
                {
                    Console.WriteLine(array.GetValue(i));
                }
                break;
            case 2:
                for (int i = 0; i < array.GetLength(0); i++)
                {
                    for (int x = 0; x < array.GetLength(1); x++)
                    {
                        Console.Write(array.GetValue(i, x));
                    }
                    Console.WriteLine();
                }
                break;
        }
    }

        private void button3_Click(object sender, EventArgs e)
        {
            // Create dictionary and add five keys and values.
            var dictionary = new Dictionary<string, int>();
            dictionary.Add("car", 2);
            dictionary.Add("apple", 1);
            dictionary.Add("zebra", 0);
            dictionary.Add("mouse", 5);
            dictionary.Add("year", 3);

            // Acquire keys and sort them.
            var list = dictionary.Keys.ToList();
            list.Sort();

            // Loop through keys.
            foreach (var key in list)
            {
                Console.WriteLine("{0}: {1}", key, dictionary[key]);
            }

            string[] colors = new string[]
        {
            "orange",
            "blue",
            "yellow",
            "aqua",
            "red"
        };
            // Call Array.Sort method.
            Array.Sort(colors);
            foreach (string color in colors)
            {
                Console.WriteLine(color);
            }
        }

        private void button4_Click(object sender, EventArgs e)
        {
            double d = Math.Round(1.3665, 1);
            Console.WriteLine(d);
            MessageBox.Show(d.ToString());

            double before1 = 123.46;
            double after1 = Math.Round(before1, 1);
            MessageBox.Show(after1.ToString());

            /*
             >>> np.around([0.37, 1.64])
array([ 0.,  2.])
>>> np.around([0.37, 1.64], decimals=1)
array([ 0.4,  1.6])
>>> np.around([.5, 1.5, 2.5, 3.5, 4.5]) # rounds to nearest even value
array([ 0.,  2.,  2.,  4.,  4.])
>>> np.around([1,2,3,11], decimals=1) # ndarray of ints is returned
array([ 1,  2,  3, 11])
>>> np.around([1,2,3,11], decimals=-1)
array([ 0,  0,  0, 10])
*/
           
        }

        private void button5_Click(object sender, EventArgs e)
        {
            /*
            var mySequence = Enumerable.Range(0, 12);
            var seq = Enumerable.Range(0, 12).ToList();
            Console.WriteLine(mySequence);
            MessageBox.Show(seq.ToString());

            IEnumerable<double> squares = (IEnumerable<double>)Enumerable.Range(1, 100);
            //Py.RangeExcl(12, 18);    // [12, 13, 14, 15, 16, 17]

            //Py.RangeIncl(12, 18);    // [12, 13, 14, 15, 16, 17, 18]
            // 
           
            IEnumerable<int> squaes = Enumerable.Range(2, 8);
        
            foreach (int num in squaes)
            {
                Console.WriteLine(num);
                MessageBox.Show((num).ToString());
            }
             */
            foreach (var index in Range(5, 7))
            {
                Console.WriteLine(index);
            }

            /*
             numpy.arange([start, ]stop, [step, ]dtype=None)
            >>> np.arange(3)
array([0, 1, 2])
>>> np.arange(3.0)
array([0., 1., 2.])
>>> np.arange(3, 7)
array([3, 4, 5, 6])
>>> np.arange(3, 7, 2)
array([3, 5])
*/
        }

        private void button6_Click(object sender, EventArgs e)
        {

            var np = new[] { 1, 2, 4, 7, 0 };
            var res = np.Zip(np.Skip(1), (a, b) => b - a).ToArray();
            foreach (var n in res)
            {
                Console.WriteLine(n);
            }

                /*
                >>> x = np.array([1, 2, 4, 7, 0])
    >>> np.diff(x)
    array([1, 2, 3, -7])
    >>> np.diff(x, n = 2)
    array([1, 1, -10])
    */
            }

        private void button7_Click(object sender, EventArgs e)
        {
            int[] array1 = { 1, 2, 3, 4 };
            int[] array2 = { 5, 6, 7, 8 };

            int[] array3 = new int[array1.Length + array2.Length];
            array1.CopyTo(array3, 0);
            array2.CopyTo(array3, array1.Length);
            Console.WriteLine(array3.Length);
              //Console.Write(i);
                Console.WriteLine(string.Join(",", array3));
            


            /*
             * 
             * int[] array1 = { 1, 3, 5 };
            int[] array2 = { 0, 2, 4 };

            // Concatenate array1 and array2.
            var result1 = array1.Concat(array2);

            Console.WriteLine(result1);


             *  int[] x = new int[] { 1, 2, 3 };
            int[] y = new int[] { 4, 5 };

            var z = new int[x.Length + y.Length];
            x.CopyTo(z, 0);
            y.CopyTo(z, x.Length);
            Console.WriteLine(z);


             >>> a = np.array((1,2,3))
>>> b = np.array((2,3,4))
>>> np.hstack((a,b))
array([1, 2, 3, 2, 3, 4])
>>> a = np.array([[1],[2],[3]])
>>> b = np.array([[2],[3],[4]])
>>> np.hstack((a,b))
array([[1, 2],
       [2, 3],
       [3, 4]])
       */
        }
    }
}
