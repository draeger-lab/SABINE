package cv;

import java.io.File;

public class DirectoryManager {

    String indent = "";

    /**
     * Prints the contents of his directory
     *
     * @param dir
     */
    public void enterDirectory(File dir) {
        System.out.println(indent + "[" + dir.getName() + "]");
        indent += "  ";
    }

    /**
     * Deletes a directory tree
     *
     * @param path
     */
    private static void deleteTree(File path) {
        if (path.listFiles() != null) {
            for (File file : path.listFiles()) {
                if (file.isDirectory())
                    deleteTree(file);
                file.delete();
            }
            path.delete();
        }
    }

    /**
     * Deletes a single directory
     *
     * @param path
     */
    public static void deleteDir(String path) {
        deleteTree(new File(path));
    }

    /**
     * returns a list f all files in this directory
     *
     * @param dir
     * @return list of files in this directory
     */
    public String[] getDirectoryContent(File dir) {
        indent += "  ";
        return dir.list();
    }

    /**
     * creates a directory
     *
     * @param path
     */
    public static void createDirectory(String path) {
        File f = new File(path);
        f.mkdirs();
    }

    /**
     * cleans a directory of all files, except the one with this name
     *
     * @param path
     * @param f
     */
    public static void cleanDirExceptThisFile(File path, String f) {

        if (path.listFiles() != null) {
            for (File file : path.listFiles()) {
                if (file.isFile() && !(file.getName().toLowerCase().equals(f.toLowerCase()))) {
                    file.delete();
                }
            }
        }
    }


    /**
     * Renames a file
     *
     * @param fIn
     * @param fOut
     */
    public static void renameThisFile(File fIn, File fOut) {
        if (fIn != null) {
            fIn.renameTo(fOut);
        }
    }
}