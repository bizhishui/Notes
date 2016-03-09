% Notes in learning Java
% Jinming LYU
% 9 mars 2016

# Java Serializable #
序列化的过程就是对象写入字节流和从字节流中读取对象。将对象状态转换成字节流之后，可以用java.io包中的各种字节流类
将其保存到文件中，管道到另一线程中或通过网络连接将对象数据发送到另一主机。对象序列化功能非常简单、强大，
在RMI、Socket、JMS、EJB都有应用。对象序列化问题在网络编程中并不是最激动人心的课题，但却相当重要，具有许多实用意义。

1. 对象序列化可以实现分布式对象。主要应用例如：RMI要利用对象序列化运行远程主机上的服务，就像在本地机上运行对象时一样。

2. java对象序列化不仅保留一个对象的数据，而且递归保存对象引用的每个对象的数据。可以将整个对象层次写入字节流中，
可以保存在文件中或在网络连接上传递。利用对象序列化可以进行对象的“深复制”，即复制对象本身及引用的对象本身。
序列化一个对象可能得到整个对象序列。

对象序列化将数据分解成字节流，以便存储在文件中或在网络上传输。反序列化就是打开字节流并重构对象。
对象序列化不仅要将基本数据类型转换成字节表示，有时还要恢复数据。恢复数据要求有恢复数据的对象实例。

## 序列化是什么？ ##
序列化就是将一个 **对象** 的状态（各个属性量）保存起来，然后在适当的时候再获得。序列化分为两大部分：序列化和反序列化。
序列化是这个过程的第一部分，将数据分解成字节流，以便存储在文件中或在网络上传输。反序列化就是打开字节流并重构对象。
对象序列化不仅要将基本数据类型转换成字节表示，有时还要恢复数据。恢复数据要求有恢复数据的对象实例。
简单说就是为了保存在内存中的各种对象的状态（也就是实例变量，不是方法），并且可以把保存的对象状态再读出来。
虽然你可以用你自己的各种各样的方法来保存object states，但是Java给你提供一种应该比你自己好的保存对象状态的机制，那就是序列化。

## 序列化的[特点](http://www.oschina.net/question/4873_23270) ##
如果某个类能够被序列化，其子类也可以被序列化。声明为static和transient类型的成员数据不能被序列化。
因为static代表 **类** 的状态， transient代表对象的临时数据。

## [什么情况下需要序列化](http://blog.csdn.net/fenglibing/article/details/8905490) ##
- 当你想把的内存中的对象状态保存到一个文件中或者数据库中时候；
- 当你想用套接字在网络上传送对象的时候；
- 当你想通过RMI传输对象的时候；

## 相关注意事项  ##
- 序列化时，只对对象的状态进行保存，而不管对象的方法；
- 当一个父类实现序列化，子类自动实现序列化，不需要显式实现Serializable接口；
- 当一个对象的实例变量引用其他对象，序列化该对象时也把引用对象进行序列化；
- 并非所有的对象都可以序列化。

## 序列化和外部化的主要区别 ##
外部化和序列化是实现同一目标的两种不同方法。通过Serializable接口对对象序列化的支持是内建于核心 API 的，但是`java.io.Externalizable`
的所有实现者必须提供读取和写出的实现。Java 已经具有了对序列化的内建支持，也就是说只要制作自己的类`java.io.Serializable`，
Java 就会试图存储和重组你的对象。如果使用外部化，你就可以选择完全由自己完成读取和写出的工作。
序列化会自动存储必要的信息，用以反序列化被存储的实例，而外部化则只保存被存储的类的标识。当你通过`java.io.Serializable`
接口序列化一个对象时，有关类的信息，比如它的属性和这些属性的类型，都与实例数据一起被存储起来。在选择走Externalizable这条路时，
Java 只存储有关每个被存储类型的非常少的信息。

## 处理对象流实例 ##
java.io包有两个序列化对象的类。ObjectOutputStream负责将对象写入字节流，ObjectInputStream从字节流重构对象。
writeObject()方法是最重要的方法，用于对象序列化。如果对象包含其他对象的引用，则writeObject()方法递归序列化这些对象。
每个ObjectOutputStream维护序列化的对象引用表，防止发送同一对象的多个拷贝。（这点很重要）由于writeObject()可以序列化整组交叉引用的对象，
因此同一ObjectOutputStream实例可能不小心被请求序列化同一对象。这时，进行反引用序列化，而不是再次写入对象字节流。

如序列化 today’s date 到一个文件中

**`FileOutputStream f = new FileOutputStream(“tmp”);`//创建一个包含恢复对象(即对象进行反序列化信息)的”tmp”数据文件**

**`ObjectOutputStream s = new ObjectOutputStream(f);`**

**`s.writeObject(“Today”);` //写入字符串对象;**

**`s.writeObject(new Date());` //写入瞬态对象;**

**`s.flush();`**

ObjectInputStream类与ObjectOutputStream相似，
它扩展DataInput接口。`ObjectInputStream`中的方法镜像DataInputStream中读取Java基本数据类型的公开方法。
`readObject()`方法从字节流中反序列化对象。每次调用`readObject()`方法都返回流中下一个Object。
对象字节流并不传输类的字节码，而是包括类名及其签名。`readObject()`收到对象时，JVM装入头中指定的类。如果找不到这个类，
则`readObject()`抛出`ClassNotFoundException`,如果需要传输对象数据和字节码，则可以用RMI框架。
`ObjectInputStream`的其余方法用于定制反序列化过程。

如从文件中反序列化 string 对象和 date 对象。

**`FileInputStream in = new FileInputStream(“tmp”);`**

**`ObjectInputStream s = new ObjectInputStream(in);`**

**`String today = (String)s.readObject();` //恢复对象;**

**`Date date = (Date)s.readObject();`**
