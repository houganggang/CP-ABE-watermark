package com.duwei.cp.abe.test;

import com.duwei.cp.abe.attribute.Attribute;
import com.duwei.cp.abe.engine.CpAneEngine;
import com.duwei.cp.abe.parameter.*;
import com.duwei.cp.abe.structure.AccessTree;
import com.duwei.cp.abe.structure.AccessTreeBuildModel;
import com.duwei.cp.abe.structure.AccessTreeNode;
import com.duwei.cp.abe.text.CipherText;
import com.duwei.cp.abe.text.PlainText;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.List;

/**
 * @BelongsProject: CP-ABE
 * @BelongsPackage: com.duwei.cp.abe.text
 * @Author: duwei
 * @Date: 2022/7/25 16:32
 * @Description: 测试类
 */
public class Test {


    public static void test1() throws IOException {
        //1.生成系统密钥  包含公钥-私钥
        SystemKey systemKey = SystemKey.build();
        //2.设置用户属性
        List<Attribute> attributes = Arrays.asList(
                new Attribute("学生", systemKey.getPublicKey()),
                new Attribute("老师", systemKey.getPublicKey()),
                new Attribute("硕士", systemKey.getPublicKey()),
                new Attribute("护士", systemKey.getPublicKey()),
                new Attribute("二班", systemKey.getPublicKey())
        );
        //3.生成用户私钥
        CpAneEngine cpAneEngine = new CpAneEngine();
        UserPrivateKey userPrivateKey = cpAneEngine.keyGen(systemKey.getMasterPrivateKey(), attributes);

        //公钥和私钥的写入
        String fileName_SystemKey = "systemKey.txt";
        String fileName_userPrivateKey = "userPrivateKey.txt";
        writeFile(fileName_SystemKey, systemKey.toString());
        writeFile(fileName_userPrivateKey, userPrivateKey.toString());

        //4.构建访问树
        AccessTree accessTree = getAccessTree(systemKey.getPublicKey());

        //5. 调用加密方法
        function(cpAneEngine, systemKey, userPrivateKey, accessTree);

//        //6.加密
//        CipherText cipherText = cpAneEngine.encrypt(systemKey.getPublicKey(), plainText, accessTree);
//        System.out.println("cipherText : " + cipherText);
//
//        //7.保存加密数据到文件
//        String fileName_cipherText = "cipherText.txt";
//        writeFile(fileName_cipherText, cipherText.toString());
//
//        //8.解密
//        String decryptStr = cpAneEngine.decryptToStr(systemKey.getPublicKey(), userPrivateKey, cipherText);
//        System.out.println("decryptStr : " + decryptStr);
//
//        //9.保存解密数据到文件
//        String fileName_decryptStrText = "decryptStrText.txt";
//        writeFile(fileName_decryptStrText, decryptStr);
    }

    public static void main(String[] args) throws IOException {
        test1();
    }


    public static AccessTree getAccessTree(PublicKey publicKey) {
        AccessTreeBuildModel[] accessTreeBuildModels = new AccessTreeBuildModel[7];
        //根节点ID必须为1
        accessTreeBuildModels[0] = AccessTreeBuildModel.innerAccessTreeBuildModel(1, 2, 1, -1);
        accessTreeBuildModels[1] = AccessTreeBuildModel.leafAccessTreeBuildModel(2, 1, "学生", 1);
        accessTreeBuildModels[2] = AccessTreeBuildModel.leafAccessTreeBuildModel(3, 2, "老师", 1);
        accessTreeBuildModels[3] = AccessTreeBuildModel.leafAccessTreeBuildModel(4, 3, "硕士", 1);
        accessTreeBuildModels[4] = AccessTreeBuildModel.innerAccessTreeBuildModel(5, 1, 4, 1);
        accessTreeBuildModels[5] = AccessTreeBuildModel.leafAccessTreeBuildModel(6, 1, "二班", 5);
        accessTreeBuildModels[6] = AccessTreeBuildModel.leafAccessTreeBuildModel(7, 2, "护士", 5);
        return AccessTree.build(publicKey, accessTreeBuildModels);
    }

    public static Pairing getPairing() {
        return PairingFactory.getPairing("params/curves/a.properties");
    }

    public static void pre(AccessTreeNode node) {
        System.out.println(node);
        for (AccessTreeNode child : node.getChildren()) {
            pre(child);
        }
    }

//    /**
//     * 读取需要加密的文件内容
//     *
//     * @param filePath
//     * @return
//     * @throws IOException
//     */
//    public static String readFile(String filePath) throws IOException {
//        StringBuilder content = new StringBuilder();
//        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
//            String line;
//            while ((line = reader.readLine()) != null) {
//                content.append(line).append("\n");
//            }
//        }
//        return content.toString().trim();
//    }

    /**
     * 写入文件方法
     *
     * @param filePath
     * @param content
     * @throws IOException
     */
    private static void writeFile(String filePath, String content) throws IOException {
        // 创建文件对象
        File file = new File(filePath);

        // 如果文件不存在，则创建文件
        if (!file.exists()) {
            file.createNewFile();
        }

        // 写入内容
        try (java.io.FileWriter writer = new java.io.FileWriter(file)) {
            writer.write(content);
        }
    }

    /**
     * 加密解密方法
     *
     * @param cpAneEngine
     * @param systemKey
     * @param userPrivateKey
     * @param accessTree
     * @throws IOException
     */
    private static void function(CpAneEngine cpAneEngine, SystemKey systemKey, UserPrivateKey userPrivateKey, AccessTree accessTree) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader("3DModel.txt"));
             BufferedWriter cipherWriter = new BufferedWriter(new FileWriter("cipherText.txt"));
             BufferedWriter decryptWriter = new BufferedWriter(new FileWriter("decryptStrText.txt"))) {

            String line;

            while ((line = reader.readLine()) != null) {
                // 读取每一行，分割为单独的浮点数
                String[] parts = line.trim().split("\\s+");
                if (parts.length == 3) {
                    StringBuilder encryptedLine = new StringBuilder();
                    StringBuilder decryptedLine = new StringBuilder();

                    for (String part : parts) {
                        try {
                            // 将浮点数加密
                            float value = Float.parseFloat(part);
                            PlainText plainText = new PlainText(Float.toString(value), systemKey.getPublicKey());
                            CipherText cipherText = cpAneEngine.encrypt(systemKey.getPublicKey(), plainText, accessTree);

                            // 将加密后的数据写入密文文件
                            encryptedLine.append(cipherText.toString()).append(" ");

                            // 写入密文文件
                            cipherWriter.write(encryptedLine.toString().trim());
                            cipherWriter.newLine();

                            // 解密
                            String decryptedValue = cpAneEngine.decryptToStr(systemKey.getPublicKey(), userPrivateKey, cipherText);
                            decryptedLine.append(decryptedValue).append(" ");

                        } catch (Exception e) {
                            System.err.println("Error processing value: " + part);
                            e.printStackTrace();
                        }
                    }

                    // 写入解密文件（保持格式）
                    decryptWriter.write(decryptedLine.toString().trim());
                    decryptWriter.newLine();
                }
            }
        } catch (IOException e) {
            System.err.println("Error reading or writing files.");
            e.printStackTrace();
        }
    }
}
