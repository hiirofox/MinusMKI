/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
LModelAudioProcessorEditor::LModelAudioProcessorEditor(LModelAudioProcessor& p)
	: AudioProcessorEditor(&p), audioProcessor(p)
{
	// Make sure that before the constructor has finished, you've set the
	// editor's size to whatever you need it to be.
	setResizable(true, true); // 允许窗口调整大小

	setOpaque(false);  // 允许在边框外面绘制

	//setResizeLimits(64 * 11, 64 * 5, 10000, 10000); // 设置最小宽高为300x200，最大宽高为800x600
	setSize(64 * 8, 64 * 2);
	setResizeLimits(64 * 8, 64 * 2, 64 * 13, 64 * 2);

	//constrainer.setFixedAspectRatio(11.0 / 4.0);  // 设置为16:9比例
	//setConstrainer(&constrainer);  // 绑定窗口的宽高限制

	K_Pitch.setText("pitch", "");
	K_Pitch.ParamLink(audioProcessor.GetParams(), "pitch");
	addAndMakeVisible(K_Pitch);
	K_P1.setText("sync", "");
	K_P1.ParamLink(audioProcessor.GetParams(), "p1");
	addAndMakeVisible(K_P1);
	K_P2.setText("pwm", "");
	K_P2.ParamLink(audioProcessor.GetParams(), "p2");
	addAndMakeVisible(K_P2);
	K_P3.setText("form", "");
	K_P3.ParamLink(audioProcessor.GetParams(), "p3");
	addAndMakeVisible(K_P3);
	K_P4.setText("fb", "");
	K_P4.ParamLink(audioProcessor.GetParams(), "p4");
	addAndMakeVisible(K_P4);
	K_P5.setText("detune", "");
	K_P5.ParamLink(audioProcessor.GetParams(), "p5");
	addAndMakeVisible(K_P5);
	K_P6.setText("param6", "");
	K_P6.ParamLink(audioProcessor.GetParams(), "p6");
	addAndMakeVisible(K_P6);


	startTimerHz(30);

}

LModelAudioProcessorEditor::~LModelAudioProcessorEditor()
{
}

//==============================================================================
void LModelAudioProcessorEditor::paint(juce::Graphics& g)
{
	g.fillAll(juce::Colour(0x00, 0x00, 0x00));

	g.fillAll();
	g.setFont(juce::Font("FIXEDSYS", 17.0, 1));
	g.setColour(juce::Colour(0xff00ff00));;

	int w = getBounds().getWidth(), h = getBounds().getHeight();

	//g.drawText("L-MODEL Magnitudelay", juce::Rectangle<float>(32, 16, w, 16), 1);
}

void LModelAudioProcessorEditor::resized()
{
	juce::Rectangle<int> bound = getBounds();
	int x = bound.getX(), y = bound.getY(), w = bound.getWidth(), h = bound.getHeight();
	auto convXY = juce::Rectangle<int>::leftTopRightBottom;

	K_Pitch.setBounds(32 + 64 * 0, 32, 64, 64);
	K_P1.setBounds(32 + 64 * 1, 32, 64, 64);
	K_P2.setBounds(32 + 64 * 2, 32, 64, 64);
	K_P3.setBounds(32 + 64 * 3, 32, 64, 64);
	K_P4.setBounds(32 + 64 * 4, 32, 64, 64);
	K_P5.setBounds(32 + 64 * 5, 32, 64, 64);
	K_P6.setBounds(32 + 64 * 6, 32, 64, 64);
}

void LModelAudioProcessorEditor::timerCallback()
{
	repaint();
}
