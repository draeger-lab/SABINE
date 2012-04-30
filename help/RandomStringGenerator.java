/*
    SABINE predicts binding specificities of transcription factors.
    Copyright (C) 2009 ZBIT, University of Tübingen, Johannes Eichner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package help;

import java.io.IOException;
import java.util.Random;

public class RandomStringGenerator {


	public static String randomString(int len) {
		
		char[] abc = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z' };
		
		Random rand_gen = new Random();
		String rand_str = "";
		
		for (int i=0; i<len; i++) {
			rand_str += abc[rand_gen.nextInt(26)];
		}
		return rand_str;
	}
	
	public static int randomNumber(int len) {
		
		Random rand_gen = new Random();
		int rand_num = 0;
		
		for (int i=0; i<len; i++) {
			rand_num *= 10;
			rand_num += rand_gen.nextInt(10);
		}
		return rand_num;
	}
	
	public static void main(String[] args) throws IOException {
		String test = randomString(10);
		System.out.println(test);
	}

}

