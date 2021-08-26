public class Test {
    public static void main(String[] args) {

        BCH bch2 = new BCH(15, 2);
        bch2.get_params_info();
        System.out.print("Generating polynomial: ");
        BCH.get_polynomial_representation(bch2.generating_polynomial);
        int[] message1 = new int[]{1, 0, 1, 0, 1, 0, 1};
        System.out.println("Input message: ");
        for (int i = 0; i < message1.length; i++) {
            System.out.print(message1[i] + " ");
        }
        System.out.println();
        int[] code1 = bch2.coding(message1);
        System.out.print("Coded input message: ");
        BCH.get_polynomial_representation(code1);
        int[] distorted_code = bch2.set_distortion(code1, new int[]{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1});
        System.out.print("Distorted code: ");
        BCH.get_polynomial_representation(distorted_code);
        System.out.println("Starting the process of decoding.");
        int[] decoded_code1 = bch2.decoding(distorted_code);
        System.out.println("Decoded message: ");
        for (int i = 0; i < decoded_code1.length; i++) {
            System.out.print(decoded_code1[i] + " ");
        }
        System.out.println('\n');

        BCH bch3 = new BCH(31, 3);
        bch3.get_params_info();
        System.out.print("Generating polynomial: ");
        BCH.get_polynomial_representation(bch3.generating_polynomial);
        int[] message2 = new int[]{1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        System.out.println("Input message: ");
        for (int i = 0; i < message2.length; i++) {
            System.out.print(message2[i] + " ");
        }
        System.out.println();

        int[] code3 = bch3.coding(message2);
        System.out.print("Coded input message: ");
        BCH.get_polynomial_representation(code3);
        int[] distortion2 = new int[]{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1};
        System.out.print("Distortion: ");
        BCH.get_polynomial_representation(distortion2);
        int[] distorted_code3 = bch3.set_distortion(code3, distortion2);
        System.out.print("Distorted code: ");
        BCH.get_polynomial_representation(distorted_code3);
        int[] decoded_message2 = bch3.decoding(distorted_code3);
        System.out.println("Decoded message:");
        for (int i = 0; i < decoded_message2.length; i++) {
            System.out.print(decoded_message2[i] + " ");
        }
        System.out.println();
    }
}
