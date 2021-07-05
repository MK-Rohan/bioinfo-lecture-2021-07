#! /usr/bin/env python

class Cat:
    def __init__(self):
        self.sleepy = True

    def snack(self):
        print('myeo~')

class KoShort(Cat):
    def snack(self):
        print("야옹")
    def setAge(self, Age):
        self.__age = Age
        print('set age to ', Age)
    def showAge(self):
        print(self.__age, 'years old.')


catcat = Cat()
catcat.snack()
print(catcat.sleepy, end = '\n\n')

minyong = KoShort()
minyong.snack()
minyong.setAge(7)
minyong.showAge()
print(minyong.sleepy)

print(minyong.__age) # 비공개 속성에 직접 접근하는 명령. >> 비공개기 때문에 에러남
