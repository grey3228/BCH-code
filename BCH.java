import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

class BCH {

    // -- parameters -- //

    private int n;  // code whole length
    private int[] possible_lengths = {3, 7, 15, 31, 63, 127, 255, 511, 1023};
    private int m;

    private int k;  // info digits length
    private int r;  // checking digits length
    private int t;  // digits to correct

    private int d;  // min code distance

    // --   --  //

    private int[][] primitive_polynomials = {{1, 1, 1}, {1, 1, 0, 1}, {1, 1, 0, 0, 1}, {1, 0, 1, 0, 0, 1}, {1, 1, 0, 0, 0, 0, 1}, {1, 0, 0, 1, 0, 0, 0, 1}, {1, 0, 1, 1, 1, 0, 0, 0, 1}, {1, 0, 0, 0, 1, 0, 0, 0, 0, 1}, {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}};    //irresolvable polynomials for constructing Galois fields from degree m = 2 up to m = 10. Fields, based on these polynomials, definitely have primitive element {x} = {0, 1}
    public int[][] primitive_element_in_degrees;
    public Integer[] generating_polynomial;

    public BCH(int n, int t) {
        // checking correctness of input parameters
        if (n > 1023) {
            System.out.println("Too long code length");
            System.exit(1);
        }
        if (n < 3) {
            System.out.println("Incorrect code length");
            System.exit(1);
        }
        boolean check_completed = false;
        for (int i = 0; i < possible_lengths.length; i++) {
            if (possible_lengths[i] == n) {

                this.m = i + 2;
                if (t >= Math.pow(2.0, (double) (m - 1))) {
                    System.out.println("Too many errors to correct");
                    System.exit(1);
                }
                this.r = this.m * t;
                this.t = t;
                this.n = n;
                this.k = this.n - this.r;
                this.d = 2 * this.t + 1;    // d >= 2 * t + 1
                check_completed = true;
                break;
            }
        }
        if (check_completed == false) {
            System.out.println("Incorrect input n");
            System.exit(1);
        }

        // We need to find roots of generating polynomial
        // To do it, we construct generating polynomial by finding cyclotomic classes of generating polinomial (d - 1) roots
        int[] primitive_element = {0, 1};
        this.primitive_element_in_degrees = new int[(int) Math.pow(2.0, (double) this.m) - 1][];
        this.primitive_element_in_degrees[0] = new int[]{1};
        this.primitive_element_in_degrees[1] = primitive_element;
        for (int i = 2; i < this.primitive_element_in_degrees.length; i++) {
            this.primitive_element_in_degrees[i] = this.get_polynomial_remainder(polynomial_multiplication(this.primitive_element_in_degrees[i - 1], primitive_element), this.primitive_polynomials[this.m - 2]);
        }
        int number_of_generating_polynomial_roots = this.d - 1;
        this.generating_polynomial = this.get_generating_polynomial(this.get_cyclotomic_classes_degrees(primitive_element, number_of_generating_polynomial_roots));
    }

    // systematic coding
    public int[] coding(int[] info_message) {
        if (info_message.length == 0 || info_message.length != this.k) {
            System.out.println("Incorrect informational block length");
            System.exit(1);
        }

        int[] output_code;
        int[] offset_polynomial = new int[this.r + 1];
        offset_polynomial[offset_polynomial.length - 1] = 1;
        output_code = this.polynomial_addition(this.polynomial_multiplication(offset_polynomial, info_message), this.get_polynomial_remainder(this.polynomial_multiplication(offset_polynomial, info_message), this.generating_polynomial));
        return output_code;
    }

    public boolean isNull(Integer[] arr) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == 1) return false;
        }
        return true;
    }

    public boolean isNull(int[] arr) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == 1) return false;
        }
        return true;
    }

    public int[] set_distortion(int[] code, int[] distortion) {
        if (code.length > this.n) {
            System.out.println("Incorrect code length");
            System.exit(1);
        }

        if (distortion.length > this.n) {
            System.out.println("Incorrect distortion length");
        }

        int error_counter = 0;
        for (int i = 0; i < distortion.length; i++) {
            if (distortion[i] == 1) error_counter++;
        }
        if (error_counter > this.t) System.out.println("Too many errors are given. This code won't fix them.");
        return this.polynomial_addition(code, distortion);
    }

    public int[] decoding(int[] code) {
        if (code.length > this.n) {
            System.out.println("Incorrect code to decode: too long code.");
            System.exit(1);
        }
        Integer[] syndroms = this.count_syndroms(code);
        for (int i = 0; i < syndroms.length; i++) {
            System.out.println("Syndrom #" + (i + 1) + ":" + syndroms[i]);
        }
        Integer[][] locators_polynomial = Berlekamp_Massey(syndroms);

        // completing Chien's search to find roots - positions of errors
        int[] error_pointers = Chien_search(locators_polynomial);
        for (int i = 0; i < error_pointers.length; i++) {
            System.out.println("Error #" + (i + 1) + " position is: " + error_pointers[i]);
        }
        // correcting errors
        int max_pointer = 0;
        for (int p = 0; p < error_pointers.length; p++) {
            if (error_pointers[p] > max_pointer) max_pointer = error_pointers[p];
        }
        int[] corrective_vector = new int[max_pointer + 1];
        for (int p = 0; p < error_pointers.length; p++) {
            corrective_vector[error_pointers[p]] = 1;
        }
        int[] corrected_code = this.polynomial_addition(corrective_vector, code);
        System.out.println("Corrected code: ");
        BCH.get_polynomial_representation(corrected_code);
        int[] info_block = new int[this.k];
        for (int p = 0; p < this.k; p++) {
            try {
                info_block[p] = corrected_code[p + this.r];
            } catch (ArrayIndexOutOfBoundsException ex) {
                info_block[p] = 0;
            }
        }
        return info_block;
    }

    public int[] Chien_search(Integer[][] locators_polynomial) {
        int[] reversed_roots = new int[get_degree(locators_polynomial[1])];
        boolean add_one;
        if (locators_polynomial[1][0] == 1) add_one = true;
        else add_one = false;
        int i = 0;
        for (int deg = 0; deg < ((int) Math.pow(2.0, (int) this.m) - 1); deg++) {

            int[] sum = new int[]{0};
            if (add_one == true) sum = new int[]{1};
            for (int x_deg = 1; x_deg < locators_polynomial[1].length; x_deg++) {

                if (locators_polynomial[1][x_deg] == 1) {
                    int curr_term_deg;
                    if (deg == 0) curr_term_deg = 0;
                    else curr_term_deg = (deg * x_deg) % ((int) Math.pow(2.0, (int) this.m) - 1);
                    curr_term_deg = (curr_term_deg + locators_polynomial[0][x_deg]) % ((int) Math.pow(2.0, (int) this.m) - 1);
                    sum = this.polynomial_addition(sum, this.primitive_element_in_degrees[curr_term_deg]);
                }
            }
            if (isNull(sum) == true) {
                for (int k = 0; k < ((int) Math.pow(2.0, (int) this.m) - 1); k++) {
                    if ((deg + k) % ((int) Math.pow(2.0, (int) this.m) - 1) == 0) {
                        reversed_roots[i] = k;
                        i++;
                        break;
                    }
                }
            }
        }
        return reversed_roots;
    }

    public Integer[][] Berlekamp_Massey(Integer[] syndroms) {
        // default values
        Map<Integer, Integer[][]> locators = new HashMap<Integer, Integer[][]>();
        Map<Integer, Integer> l = new HashMap<Integer, Integer>();  // l_{-1}, l_{0}
        l.put((-1), 0);
        l.put(0, 0);
        Map<Integer, Integer> d = new HashMap<Integer, Integer>();    // in fact, it stores degrees of primitive element
        d.put((-1), 0);   // epsilon^{0} = 1
        if (syndroms[0] != (-1)) d.put(0, syndroms[0]);
        else d.put(0, -1);
        locators.put((-1), new Integer[][]{{0}, {1}});
        locators.put(0, new Integer[][]{{0}, {1}});
        int i;
        for (i = 0; ; i++) {
            // case 1
            if (d.get(i) == -1) { // means d = 0
                locators.put(i + 1, locators.get(i));
                l.put(i + 1, l.get(i));

                if (i >= (l.get(i + 1) + this.t - 1) || (i == 2 * this.t - 1)) break;
                int d_i_next;
                d_i_next = syndroms[i + 1];

                for (int k = i, u = 1; k >= i - l.get(i + 1) + 1 && u <= l.get(i + 1); k = k - 1, u = u + 1) {
                    int curr_addition;
                    if (syndroms[k] != -1 && locators.get(i + 1)[0][u] != -1) {
                        curr_addition = (syndroms[k] + locators.get(i + 1)[0][u]) % ((int) Math.pow(2.0, (int) this.m) - 1);
                    } else curr_addition = -1;
                    int[] curr_sum;
                    if (d_i_next != -1 && curr_addition != -1) {
                        curr_sum = this.polynomial_addition(primitive_element_in_degrees[d_i_next], this.primitive_element_in_degrees[curr_addition]);

                        int nulls = 0;
                        for (int t = curr_sum.length - 1; t >= 0; t--) {
                            if (curr_sum[t] == 1) break;
                            if (curr_sum[t] == 0) {
                                nulls++;
                            }
                        }
                        int[] curr_sum_normalized;
                        if (nulls != 0) {
                            curr_sum_normalized = new int[curr_sum.length - nulls];
                            for (int t = 0; t < curr_sum_normalized.length; t++) {
                                curr_sum_normalized[t] = curr_sum[t];
                            }
                        } else curr_sum_normalized = curr_sum;
                        d_i_next = -1;
                        if (isNull(curr_sum_normalized) != true) {
                            for (int j = 0; j < this.primitive_element_in_degrees.length; j++) {

                                if (curr_sum_normalized.length == this.primitive_element_in_degrees[j].length) {
                                    boolean found = true;
                                    for (int y = 0; y < this.primitive_element_in_degrees[j].length; y++) {
                                        if (curr_sum[y] != this.primitive_element_in_degrees[j][y]) {
                                            found = false;
                                            break;
                                        }
                                    }
                                    if (found == true) {
                                        d_i_next = j;
                                        break;
                                    }

                                }
                            }
                        }
                    } else if (d_i_next == -1 && curr_addition != -1) {
                        d_i_next = curr_addition;
                    } else if (d_i_next == -1 && curr_addition == -1) {
                        d_i_next = -1;
                    }
                }
                d.put(i + 1, d_i_next);

            } else {
                // first of all, counting m
                int m = -100;   // value by default
                int difference = -10;
                for (int m_trial = -1; m_trial >= (-1) && m_trial < i; m_trial++) {
                    if (m_trial - l.get(m_trial) > difference && (d.get(m_trial) != -1)) {
                        difference = m_trial - l.get(m_trial);
                        m = m_trial;
                    }
                }
                Integer[][] first_term = locators.get(i);
                Integer d_i = d.get(i);
                Integer d_m_reversed = 0;
                for (int p = 0; p < ((int) Math.pow(2.0, (int) this.m) - 1); p++) {
                    if ((p + d.get(m)) % ((int) Math.pow(2.0, (int) this.m) - 1) == 0) {
                        d_m_reversed = p;
                        break;
                    }
                }

                Integer[] x_i_m_diff = new Integer[i - m + 1];
                x_i_m_diff[x_i_m_diff.length - 1] = 1;
                int x_locator_mul_deg = 0;

                if (get_degree(locators.get(m)[1]) > 0 && get_degree(x_i_m_diff) > 0)
                    x_locator_mul_deg = get_degree(locators.get(m)[1]) + get_degree(x_i_m_diff);
                else if (get_degree(locators.get(m)[1]) == 0 && get_degree(x_i_m_diff) > 0)
                    x_locator_mul_deg = get_degree(x_i_m_diff);
                else if (get_degree(locators.get(m)[1]) > 0 && get_degree(x_i_m_diff) == 0)
                    x_locator_mul_deg = get_degree(locators.get(m)[1]);
                else x_locator_mul_deg = 0;
                Integer[][] second_term = new Integer[2][x_locator_mul_deg + 1];
                int offset = x_i_m_diff.length - 1;

                second_term[1] = this.polynomial_multiplication(locators.get(m)[1], x_i_m_diff);

                for (int u = 0; u < offset; u++) {
                    second_term[0][u] = -1;
                }
                for (int u = offset; u < second_term[0].length; u++) {
                    if (locators.get(m)[1][u - offset] == 1) {
                        second_term[0][u] = locators.get(m)[0][u - offset];
                    }
                }

                for (int u = 0; u < second_term[0].length; u++) {
                    if (second_term[0][u] != -1) {
                        second_term[0][u] = (second_term[0][u] + d_m_reversed + d_i) % ((int) Math.pow(2.0, (int) this.m) - 1);
                    }
                }
                // finally, summing terms:
                Integer[][] new_locator = new Integer[2][];
                new_locator = this.sum_locators(first_term, second_term);
                locators.put(i + 1, new_locator);
                l.put(i + 1, Math.max(l.get(i), l.get(m) + i - m));
                if (i >= (l.get(i + 1) + this.t - 1) || (i == 2 * this.t - 1)) break;
                int d_i_next;
                d_i_next = syndroms[i + 1];
                for (int k = i, u = 1; k >= i - l.get(i + 1) + 1 && u <= l.get(i + 1); k = k - 1, u = u + 1) {
                    int curr_addition;
                    if (syndroms[k] != -1 && locators.get(i + 1)[0][u] != -1) {
                        curr_addition = (syndroms[k] + locators.get(i + 1)[0][u]) % ((int) Math.pow(2.0, (int) this.m) - 1);
                    } else curr_addition = -1;
                    int[] curr_sum;
                    if (d_i_next != -1 && curr_addition != -1) {
                        curr_sum = this.polynomial_addition(primitive_element_in_degrees[d_i_next], this.primitive_element_in_degrees[curr_addition]);

                        int nulls = 0;
                        for (int t = curr_sum.length - 1; t >= 0; t--) {
                            if (curr_sum[t] == 1) break;
                            if (curr_sum[t] == 0) {
                                nulls++;
                            }
                        }
                        int[] curr_sum_normalized;
                        if (nulls != 0) {
                            curr_sum_normalized = new int[curr_sum.length - nulls];
                            for (int t = 0; t < curr_sum_normalized.length; t++) {
                                curr_sum_normalized[t] = curr_sum[t];
                            }
                        } else curr_sum_normalized = curr_sum;
                        d_i_next = -1;
                        if (isNull(curr_sum_normalized) != true) {
                            for (int j = 0; j < this.primitive_element_in_degrees.length; j++) {

                                if (curr_sum_normalized.length == this.primitive_element_in_degrees[j].length) {
                                    boolean found = true;
                                    for (int y = 0; y < this.primitive_element_in_degrees[j].length; y++) {
                                        if (curr_sum[y] != this.primitive_element_in_degrees[j][y]) {
                                            found = false;
                                            break;
                                        }
                                    }
                                    if (found == true) {
                                        d_i_next = j;
                                        break;
                                    }

                                }
                            }
                        }
                    } else if (d_i_next == -1 && curr_addition != -1) {
                        d_i_next = curr_addition;
                    } else if (d_i_next == -1 && curr_addition == -1) {
                        d_i_next = -1;
                    }


                }
                d.put(i + 1, d_i_next);

            }

        }
        return locators.get(i + 1);
    }

    public Integer[][] sum_locators(Integer[][] l1, Integer[][] l2) {
        // combination {-1, 0} means 0 at position j
        int new_locator_degree = Math.max(get_degree(l1[1]), get_degree(l2[1]));
        Integer[][] new_locator = new Integer[2][new_locator_degree + 1];

        int i;
        int less_length = Math.min(l1[0].length, l2[0].length);
        for (i = 0; i < less_length; i++) {
            if (l1[0][i] == -1 && l2[0][i] != -1) {
                new_locator[0][i] = l2[0][i];
                new_locator[1][i] = l2[1][i];
            } else if (l1[0][i] != -1 && l2[0][i] == -1) {
                new_locator[0][i] = l1[0][i];
                new_locator[1][i] = l1[1][i];
            } else if (l1[0][i] == -1 && l2[0][i] == -1) {
                new_locator[0][i] = -1;
                new_locator[1][i] = 0;
            } else {
                int[] new_i_eps = this.polynomial_addition(primitive_element_in_degrees[l1[0][i]], this.primitive_element_in_degrees[l2[0][i]]);
                int nulls = 0;
                for (int t = new_i_eps.length - 1; t >= 0; t--) {
                    if (new_i_eps[t] == 1) break;
                    if (new_i_eps[t] == 0) {
                        nulls++;
                    }
                }
                int[] new_i_eps_normalized;
                if (nulls != 0) {
                    new_i_eps_normalized = new int[new_i_eps.length - nulls];
                    for (int t = 0; t < new_i_eps_normalized.length; t++) {
                        new_i_eps_normalized[t] = new_i_eps[t];
                    }
                } else new_i_eps_normalized = new_i_eps;
                int new_i_eps_deg = -1;
                if (isNull(new_i_eps_normalized) != true) {
                    for (int j = 0; j < this.primitive_element_in_degrees.length; j++) {

                        if (new_i_eps_normalized.length == this.primitive_element_in_degrees[j].length) {
                            boolean found = true;
                            for (int k = 0; k < this.primitive_element_in_degrees[j].length; k++) {
                                if (new_i_eps[k] != this.primitive_element_in_degrees[j][k]) {
                                    found = false;
                                    break;
                                }
                            }
                            if (found == true) {
                                new_i_eps_deg = j;
                                break;
                            }

                        }
                    }
                }


                if (new_i_eps_deg == (-1)) {
                    new_locator[0][i] = -1;
                    new_locator[1][i] = 0;
                    continue;
                }
                new_locator[0][i] = new_i_eps_deg;
                new_locator[1][i] = 1;
            }
        }
        if (l1[0].length >= l2[0].length) {
            for (; i < l1[0].length; i++) {
                new_locator[0][i] = l1[0][i];
                new_locator[1][i] = l1[1][i];
            }
        } else {
            for (; i < l2[0].length; i++) {
                new_locator[0][i] = l2[0][i];
                new_locator[1][i] = l2[1][i];
            }
        }
        int nulls_counter = 0;
        for (int k = new_locator[1].length - 1; k >= 0; k--) {
            if (new_locator[1][k] == 1) break;
            if (new_locator[1][k] == null || new_locator[0][k] == -1) {
                nulls_counter++;
            }
        }
        if (nulls_counter == new_locator[1].length) return new Integer[][]{{-1}, {0}};
        Integer[][] normalized_locator = new Integer[2][new_locator[1].length - nulls_counter];
        for (int k = 0; k < normalized_locator[0].length; k++) {
            normalized_locator[0][k] = new_locator[0][k];
            normalized_locator[1][k] = new_locator[1][k];
        }
        return normalized_locator;
    }


    private Integer[] count_syndroms(int[] code) {
        Integer[] syndroms = new Integer[this.t * 2];
        int[][] syndroms_polynomials = new int[this.t * 2][];
        for (int i = 0; i < syndroms.length; i++) {
            syndroms_polynomials[i] = new int[]{0};
            if (code[0] == 1) syndroms_polynomials[i][0] = 1;
            for (int j = 1; j < code.length; j++) {
                if (code[j] == 1) {
                    int curr_overall_degree = ((i + 1) * (j)) % ((int) Math.pow(2.0, (int) this.m) - 1);
                    syndroms_polynomials[i] = this.polynomial_addition(syndroms_polynomials[i], this.primitive_element_in_degrees[curr_overall_degree]);
                }
            }
            // determine final syndrom degree
            int null_counter = 0;
            for (int p = syndroms_polynomials[i].length - 1; p >= 0; p--) {
                if (syndroms_polynomials[i][p] == 1) break;
                else null_counter++;
            }
            if (null_counter == syndroms_polynomials[i].length) { // it means syndrom_i = 0;
                syndroms[i] = -1;   // sign of null syndrom
            } else {
                int[] curr_normalized_pol = new int[syndroms_polynomials[i].length - null_counter];
                for (int p = 0; p < curr_normalized_pol.length; p++) {
                    curr_normalized_pol[p] = syndroms_polynomials[i][p];
                }
                for (int j = 0; j < ((int) Math.pow(2.0, (int) this.m) - 1); j++) {
                    boolean found = false;
                    if (curr_normalized_pol.length == this.primitive_element_in_degrees[j].length) {
                        found = true;
                        for (int k = 0; k < this.primitive_element_in_degrees[j].length; k++) {

                            if (curr_normalized_pol[k] != this.primitive_element_in_degrees[j][k]) {
                                found = false;
                                break;
                            }
                        }
                    }
                    if (found == true) {
                        syndroms[i] = j;
                        break;
                    }
                }
            }
        }
        return syndroms;
    }

    public int get_degree(Integer[] polynomial) {
        for (int i = polynomial.length - 1; i >= 0; i--) {
            if (polynomial[i] == 1) {
                return i;
            }
        }
        return 0;
    }

    public int get_degree(int[] polynomial) {
        for (int i = polynomial.length - 1; i >= 0; i--) {
            if (polynomial[i] == 1) {
                return i;
            }
        }
        return 0;
    }

    public Integer[] get_generating_polynomial(ArrayList<ArrayList<Integer>> cyclotomic_classes_degrees) {
        Integer[][] minimal_polynomials = new Integer[cyclotomic_classes_degrees.size()][];
        minimal_polynomials[0] = new Integer[this.primitive_polynomials[this.m - 2].length];
        for (int i = 0; i < this.primitive_polynomials[this.m - 2].length; i++) {
            minimal_polynomials[0][i] = this.primitive_polynomials[this.m - 2][i];                    // primitive polynomial has root equal primitive element, so whole cyclotomic class with primitive element represents primitive polynomial
        }

        for (int class_num = 1; class_num < cyclotomic_classes_degrees.size(); class_num++) {
            ArrayList<Integer[][]> multipliers_in_pairs = new ArrayList<Integer[][]>();  // this is a store, where we will accumulate our multipliers to sum them in the end
            // first multiplier in pair represents "x^{n}", second - "alpha^{k}"
            // multipliers in first brackets : (x + alpha^{k})
            multipliers_in_pairs.add(new Integer[2][]);
            multipliers_in_pairs.get(0)[0] = new Integer[]{0, 1};
            multipliers_in_pairs.get(0)[1] = new Integer[(int) Math.pow(2.0, (double) this.m) - 1]; // CAREFULLY!!
            // imitializing array with nulls, because Integer is a class, so by default it is a null-pointers
            for (int i = 0; i < multipliers_in_pairs.get(0)[1].length; i++) {
                multipliers_in_pairs.get(0)[1][i] = 0;
            }
            int alpha_deg_counter = 0;  // position of current degree in cyclotomic class
            multipliers_in_pairs.add(new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1]);
            for (int i = 0; i < multipliers_in_pairs.get(1)[0].length; i++) {
                multipliers_in_pairs.get(1)[0][i] = 0;
                multipliers_in_pairs.get(1)[1][i] = 0;
            }
            multipliers_in_pairs.get(1)[1][cyclotomic_classes_degrees.get(class_num).get(alpha_deg_counter)] = 1;

            for (alpha_deg_counter = 1; alpha_deg_counter < cyclotomic_classes_degrees.get(class_num).size(); alpha_deg_counter++) {
                int[][] current_bracket = new int[2][]; // (x + alpha^{cyclotomic_degree})
                current_bracket[0] = new int[]{0, 1};   // x
                current_bracket[1] = new int[cyclotomic_classes_degrees.get(class_num).get(alpha_deg_counter) + 1];
                current_bracket[1][cyclotomic_classes_degrees.get(class_num).get(alpha_deg_counter)] = 1;   // alpha ^ {cyclotomic_degree}
                int previous_multipliers_counter;

                int list_curr_len = multipliers_in_pairs.size();
                for (previous_multipliers_counter = list_curr_len - (int) Math.pow(2.0, (double) alpha_deg_counter); list_curr_len - previous_multipliers_counter > 0; previous_multipliers_counter++) {
                    int[][] new_mult_1 = new int[2][(int) Math.pow(2.0, (double) this.m) - 1], new_mult_2 = new int[2][(int) Math.pow(2.0, (double) this.m) - 1];

                    // case: (x*alpha) * (x + alpha^k)
                    if (isNull(multipliers_in_pairs.get(previous_multipliers_counter)[0]) == false && isNull(multipliers_in_pairs.get(previous_multipliers_counter)[1]) == false) {

                        new_mult_1[0] = this.polynomial_multiplication(multipliers_in_pairs.get(previous_multipliers_counter)[0], current_bracket[0]);
                        for (int i = 0; i < multipliers_in_pairs.get(previous_multipliers_counter)[1].length; i++) {
                            new_mult_1[1][i] = multipliers_in_pairs.get(previous_multipliers_counter)[1][i];
                        }

                        for (int i = 0; i < multipliers_in_pairs.get(previous_multipliers_counter)[0].length; i++) {
                            new_mult_2[0][i] = multipliers_in_pairs.get(previous_multipliers_counter)[0][i];
                        }

                        int curr_alpha_degree = 0;
                        for (int i = multipliers_in_pairs.get(previous_multipliers_counter)[1].length - 1; i > 0; i--) {
                            if (multipliers_in_pairs.get(previous_multipliers_counter)[1][i] == 1) {
                                curr_alpha_degree = i;
                                break;
                            }
                        }
                        int current_bracket_1_degree = 0;
                        for (int i = current_bracket[1].length - 1; i > 0; i--) {
                            if (current_bracket[1][i] == 1) {
                                current_bracket_1_degree = i;
                                break;
                            }
                        }
                        new_mult_2[1] = new int[((curr_alpha_degree + current_bracket_1_degree) % ((int) Math.pow(2.0, (double) this.m) - 1)) + 1];
                        new_mult_2[1][new_mult_2[1].length - 1] = 1;

                        Integer[][] cast1 = new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1];
                        Integer[][] cast2 = new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1];
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < cast1[i].length; j++) {
                                cast1[i][j] = cast2[i][j] = 0;
                            }
                        }

                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < new_mult_1[i].length; j++) {
                                cast1[i][j] = (Integer) new_mult_1[i][j];
                            }
                        }
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < new_mult_2[i].length; j++) {
                                cast2[i][j] = (Integer) new_mult_2[i][j];
                            }
                        }

                        multipliers_in_pairs.add(cast1);
                        multipliers_in_pairs.add(cast2);
                    }
                    // case: (x) * (x + alpha ^(k))
                    else if (isNull(multipliers_in_pairs.get(previous_multipliers_counter)[0]) == false && isNull(multipliers_in_pairs.get(previous_multipliers_counter)[1]) == true) {

                        new_mult_1[0] = this.polynomial_multiplication(multipliers_in_pairs.get(previous_multipliers_counter)[0], current_bracket[0]);

                        for (int i = 0; i < multipliers_in_pairs.get(previous_multipliers_counter)[0].length; i++) {
                            new_mult_2[0][i] = multipliers_in_pairs.get(previous_multipliers_counter)[0][i];

                        }

                        for (int i = 0; i < current_bracket[1].length; i++) {
                            new_mult_2[1][i] = current_bracket[1][i];
                        }

                        Integer[][] cast1 = new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1];
                        Integer[][] cast2 = new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1];
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < cast1[i].length; j++) {
                                cast1[i][j] = cast2[i][j] = 0;
                            }
                        }

                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < new_mult_1[i].length; j++) {
                                cast1[i][j] = (Integer) new_mult_1[i][j];
                            }
                        }
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < new_mult_2[i].length; j++) {
                                cast2[i][j] = (Integer) new_mult_2[i][j];
                            }
                        }

                        multipliers_in_pairs.add(cast1);
                        multipliers_in_pairs.add(cast2);

                    }

                    // case: (alpha) * (x + alpha^k)
                    else if (isNull(multipliers_in_pairs.get(previous_multipliers_counter)[0]) == true && isNull(multipliers_in_pairs.get(previous_multipliers_counter)[1]) == false) {

                        int curr_alpha_degree = 0;
                        for (int i = multipliers_in_pairs.get(previous_multipliers_counter)[1].length - 1; i > 0; i--) {
                            if (multipliers_in_pairs.get(previous_multipliers_counter)[1][i] == 1) {
                                curr_alpha_degree = i;
                                break;
                            }
                        }
                        int current_bracket_1_degree = 0;
                        for (int i = current_bracket[1].length - 1; i > 0; i--) {
                            if (current_bracket[1][i] == 1) {
                                current_bracket_1_degree = i;
                                break;
                            }
                        }
                        new_mult_2[1] = new int[((curr_alpha_degree + current_bracket_1_degree) % ((int) Math.pow(2.0, (double) this.m) - 1)) + 1];
                        new_mult_2[1][new_mult_2[1].length - 1] = 1;

                        for (int i = 0; i < current_bracket[0].length; i++) {
                            new_mult_1[0][i] = current_bracket[0][i];
                        }

                        for (int i = 0; i < multipliers_in_pairs.get(previous_multipliers_counter)[1].length; i++) {
                            new_mult_1[1][i] = multipliers_in_pairs.get(previous_multipliers_counter)[1][i];
                        }

                        Integer[][] cast1 = new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1];
                        Integer[][] cast2 = new Integer[2][(int) Math.pow(2.0, (double) this.m) - 1];
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < cast1[i].length; j++) {
                                cast1[i][j] = cast2[i][j] = 0;
                            }
                        }

                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < new_mult_1[i].length; j++) {
                                cast1[i][j] = (Integer) new_mult_1[i][j];
                            }
                        }
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < new_mult_2[i].length; j++) {
                                cast2[i][j] = (Integer) new_mult_2[i][j];
                            }
                        }

                        multipliers_in_pairs.add(cast1);
                        multipliers_in_pairs.add(cast2);
                    }
                }
            }
            int num_of_multiplications = (int) Math.pow(2.0, (double) (alpha_deg_counter));
            // getting polynomial representation of each alpha
            for (int i = multipliers_in_pairs.size() - num_of_multiplications; i < multipliers_in_pairs.size(); i++) {
                int alpha_degree = 0;
                for (int j = multipliers_in_pairs.get(i)[1].length - 1; j >= 0; j--) {
                    if (multipliers_in_pairs.get(i)[1][j] == 1) {
                        alpha_degree = j;
                        break;
                    }
                }
                alpha_degree = alpha_degree % ((int) Math.pow(2.0, (double) this.m) - 1);
                multipliers_in_pairs.get(i)[1] = new Integer[this.primitive_element_in_degrees[alpha_degree].length];
                for (int j = 0; j < this.primitive_element_in_degrees[alpha_degree].length; j++) {
                    multipliers_in_pairs.get(i)[1][j] = this.primitive_element_in_degrees[alpha_degree][j];
                }

            }

            // finally making every multiplication (x * alpha^k), where alpha^k is already represented as a polynomial
            Integer[][] polynomials_to_sum = new Integer[num_of_multiplications][];
            int counter = 0;
            for (int i = multipliers_in_pairs.size() - num_of_multiplications; i < multipliers_in_pairs.size(); i++) {
                if (isNull(multipliers_in_pairs.get(i)[0]) == false && isNull(multipliers_in_pairs.get(i)[1]) == false) {
                    polynomials_to_sum[counter] = this.polynomial_multiplication(multipliers_in_pairs.get(i)[0], multipliers_in_pairs.get(i)[1]);
                } else if (isNull(multipliers_in_pairs.get(i)[0]) == true && isNull(multipliers_in_pairs.get(i)[1]) == false) {
                    polynomials_to_sum[counter] = new Integer[multipliers_in_pairs.get(i)[1].length];
                    for (int j = 0; j < multipliers_in_pairs.get(i)[1].length; j++) {
                        polynomials_to_sum[counter][j] = multipliers_in_pairs.get(i)[1][j];
                    }
                }
                counter++;
            }

            Integer[] min_polynomial = polynomials_to_sum[0];
            for (int i = 1; i < polynomials_to_sum.length; i++) {
                min_polynomial = this.polynomial_addition(min_polynomial, polynomials_to_sum[i]);
            }

            minimal_polynomials[class_num] = new Integer[min_polynomial.length];
            for (int i = 0; i < min_polynomial.length; i++) {
                minimal_polynomials[class_num][i] = min_polynomial[i];
            }
        }

        Integer[] generating_pol = minimal_polynomials[0];
        for (int i = 1; i < cyclotomic_classes_degrees.size(); i++) {
            generating_pol = this.polynomial_multiplication(generating_pol, minimal_polynomials[i]);
        }
        return generating_pol;
    }

    public ArrayList<ArrayList<Integer>> get_cyclotomic_classes_degrees(int[] primitive_element, int highest_root_degree) {
        ArrayList<ArrayList<Integer>> cyclotomic_classes_degrees = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> collected_root_degrees = new ArrayList<Integer>();
        int p = 2;  // binary
        int mod = ((int) Math.pow(2.0, this.m)) - 1;
        int cyclotomic_classes_counter = -1;
        boolean new_cyclotomic_class_found = false;
        for (int s = 1; s <= highest_root_degree; s++) {
            int p_degree = 0;
            int curr_root_degree = (s * (int) Math.pow((double) p, (double) p_degree)) % mod;
            while (collected_root_degrees.contains(curr_root_degree) == false) {
                collected_root_degrees.add(curr_root_degree);
                if (new_cyclotomic_class_found == false) {
                    new_cyclotomic_class_found = true;
                    cyclotomic_classes_counter++;
                    cyclotomic_classes_degrees.add(new ArrayList<Integer>());
                }
                cyclotomic_classes_degrees.get(cyclotomic_classes_counter).add(curr_root_degree);
                p_degree++;
                curr_root_degree = (s * (int) Math.pow((double) p, (double) p_degree)) % mod;
            }
            new_cyclotomic_class_found = false;
        }
        return cyclotomic_classes_degrees;
    }

    public void get_params_info() {
        System.out.println("(" + this.n + ", " + this.k + ", " + this.d + ")\nn: " + this.n + "\nk: " + this.k + "\nd: " + this.d + "\nt: " + this.t + "\nm: " + this.m + "\nr: " + this.r);
    }

    public static void get_polynomial_representation(Integer[] vector) {
        int[] v = new int[vector.length];
        for (int i = 0; i < v.length; i++) {
            v[i] = vector[i];
        }
        BCH.get_polynomial_representation(v);
    }

    public static void get_polynomial_representation(int[] vector) {
        StringBuffer buff = new StringBuffer();
        int null_counter = 0;
        for (int i = 0; i < vector.length; i++) {
            if (vector[i] != 0) {
                if (i == 0) {
                    buff.append("1");
                    boolean has_next = false;
                    for (int j = 1; j < vector.length; j++) {
                        if (vector[j] != 0) {
                            has_next = true;
                            break;
                        }
                    }
                    if (has_next == true) buff.append(" + ");
                    continue;
                }
                boolean has_next = false;
                for (int j = i + 1; j < vector.length; j++) {
                    if (vector[j] != 0) {
                        has_next = true;
                        break;
                    }
                }
                if (has_next == true) buff.append("x^{" + i + "} + ");
                else buff.append("x^{" + i + "}");

            } else null_counter++;
            if (null_counter == vector.length) {
                System.out.println("0");
                return;
            }
        }
        System.out.println(buff.toString());
    }

    public int[] polynomial_addition(int[] p1, int[] p2) {
        if (p1.length >= p2.length) {
            int[] result = new int[p1.length];

            int i;
            for (i = 0; i < p2.length; i++) {
                result[i] = (p1[i] + p2[i]) % 2;
            }
            for (; i < p1.length; i++) {
                result[i] = p1[i];
            }
            return result;

        } else {
            int[] result = new int[p2.length];
            int i;
            for (i = 0; i < p1.length; i++) {
                result[i] = (p1[i] + p2[i]) % 2;
            }
            for (; i < p2.length; i++) {
                result[i] = p2[i];
            }
            return result;
        }
    }

    public Integer[] polynomial_addition(Integer[] p1, Integer[] p2) {
        if (p1.length >= p2.length) {
            Integer[] result = new Integer[p1.length];

            int i;
            for (i = 0; i < p2.length; i++) {
                result[i] = (p1[i] + p2[i]) % 2;
            }
            for (; i < p1.length; i++) {
                result[i] = p1[i];
            }
            return result;

        } else {
            Integer[] result = new Integer[p2.length];
            int i;
            for (i = 0; i < p1.length; i++) {
                result[i] = (p1[i] + p2[i]) % 2;
            }
            for (; i < p2.length; i++) {
                result[i] = p2[i];
            }
            return result;
        }
    }

    public int[] polynomial_multiplication(Integer[] p1, int[] p2) {
        int[] p_1 = new int[p1.length];
        for (int i = 0; i < p_1.length; i++) {
            p_1[i] = (int) p1[i];
        }
        return this.polynomial_multiplication(p_1, p2);
    }

    public int[] polynomial_multiplication(int[] p1, int[] p2) {
        int p1_degree = 0, p2_degree = 0;
        for (int i = p1.length - 1; i >= 0; i--) {
            if (p1[i] == 1) {
                p1_degree = i;
                break;
            }
        }
        for (int i = p2.length - 1; i >= 0; i--) {
            if (p2[i] == 1) {
                p2_degree = i;
                break;
            }
        }

        int result_length = 1;
        if (p1_degree > 0 && p2_degree > 0) result_length = p1_degree + p2_degree + 1;
        else if (p1_degree > 0 && p2_degree == 0) result_length = p1_degree + 1;
        else if (p1_degree == 0 && p2_degree > 0) result_length = p2_degree + 1;

        int[] result = new int[result_length];
        for (int i = 0; i < p1.length; i++) {
            if (p1[i] != 0) {
                int m1_degree = i;
                int[] temp_sum;
                int temp_sum_length = 1;
                if (m1_degree > 0 && p2_degree == 0) temp_sum_length = m1_degree + 1;
                if (m1_degree == 0 && p2_degree > 0) temp_sum_length = p2_degree + 1;
                if (m1_degree > 0 && p2_degree > 0) temp_sum_length = m1_degree + p2_degree + 1;
                temp_sum = new int[temp_sum_length];

                for (int j = 0; j < p2.length; j++) {
                    int m2_degree = j;
                    int multiplication_degree = m1_degree;
                    if (p2[j] != 0) {
                        multiplication_degree += m2_degree;
                        temp_sum[multiplication_degree] = 1;
                    }
                }
                result = this.polynomial_addition(result, temp_sum);
            }
        }
        return result;
    }

    public Integer[] polynomial_multiplication(Integer[] p1, Integer[] p2) {
        int p1_degree = 0, p2_degree = 0;
        for (int i = 0; i < p1.length; i++) {
            if (p1[i] == null) p1[i] = 0;
        }
        for (int i = 0; i < p2.length; i++) {
            if (p2[i] == null) p2[i] = 0;
        }
        for (int i = p1.length - 1; i >= 0; i--) {
            if (p1[i] == 1) {
                p1_degree = i;
                break;
            }
        }
        for (int i = p2.length - 1; i >= 0; i--) {
            if (p2[i] == 1) {
                p2_degree = i;
                break;
            }
        }

        int result_length = 1;
        if (p1_degree > 0 && p2_degree > 0) result_length = p1_degree + p2_degree + 1;
        else if (p1_degree > 0 && p2_degree == 0) result_length = p1_degree + 1;
        else if (p1_degree == 0 && p2_degree > 0) result_length = p2_degree + 1;

        Integer[] result = new Integer[result_length];
        for (int l = 0; l < result_length; l++) {
            result[l] = 0;
        }
        for (int i = 0; i < p1.length; i++) {
            if (p1[i] != 0) {
                int m1_degree = i;
                Integer[] temp_sum;
                int temp_sum_length = 1;
                if (m1_degree > 0 && p2_degree == 0) temp_sum_length = m1_degree + 1;
                if (m1_degree == 0 && p2_degree > 0) temp_sum_length = p2_degree + 1;
                if (m1_degree > 0 && p2_degree > 0) temp_sum_length = m1_degree + p2_degree + 1;
                temp_sum = new Integer[temp_sum_length];
                for (int l = 0; l < temp_sum_length; l++) {
                    temp_sum[l] = 0;
                }

                for (int j = 0; j < p2.length; j++) {
                    int m2_degree = j;
                    int multiplication_degree = m1_degree;
                    if (p2[j] != 0) {
                        multiplication_degree += m2_degree;
                        temp_sum[multiplication_degree] = 1;
                    }
                }
                result = this.polynomial_addition(result, temp_sum);
            }
        }
        return result;
    }

    public int[] get_polynomial_remainder(int[] p, Integer[] mod) {
        int[] mod_ = new int[mod.length];
        for (int i = 0; i < mod.length; i++) {
            mod_[i] = mod[i];
        }
        return this.get_polynomial_remainder(p, mod_);
    }

    public int[] get_polynomial_remainder(int[] p, int[] mod) {
        int p_degree = 0, mod_degree = 0;
        for (int i = p.length - 1; i >= 0; i--) {
            if (p[i] == 1) {
                p_degree = i;
                break;
            }
        }
        for (int i = mod.length - 1; i >= 0; i--) {
            if (mod[i] == 1) {
                mod_degree = i;
                break;
            }
        }
        if (p_degree < mod_degree) return p;

        int[] current_dividend = new int[p.length];
        for (int i = 0; i < p.length; i++) {
            current_dividend[i] = p[i];
        }
        int current_dividend_degree = p_degree;
        while (current_dividend_degree >= mod_degree) {
            int curr_degrees_difference = current_dividend_degree - mod_degree;
            int[] curr_multiplier = new int[curr_degrees_difference + 1];
            curr_multiplier[curr_multiplier.length - 1] = 1;
            int[] multiplied_mod = this.polynomial_multiplication(mod, curr_multiplier);
            current_dividend = this.polynomial_addition(current_dividend, multiplied_mod);

            // recount dividend degree
            current_dividend_degree = 0;
            for (int k = current_dividend.length - 1; k >= 0; k--) {
                if (current_dividend[k] == 1) {
                    current_dividend_degree = k;
                    break;
                }
            }
        }

        // slicing "0" digits
        int updated_degree = 0;
        for (int i = current_dividend.length - 1; i >= 0; i--) {
            if (current_dividend[i] == 1) {
                updated_degree = i;
                break;
            }
        }
        int[] normalized_polynomial = new int[updated_degree + 1];
        for (int i = 0; i < normalized_polynomial.length; i++) {
            normalized_polynomial[i] = current_dividend[i];
        }
        return normalized_polynomial;
    }
}